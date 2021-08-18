#include "nuchic/EventGen.hh"
#include "nuchic/Event.hh"
#include "nuchic/EventWriter.hh"
#include "nuchic/HardScatteringFactory.hh"
#include "nuchic/HardScattering.hh"
#include "nuchic/Nucleus.hh"
#include "nuchic/Beams.hh"
#include "nuchic/Cascade.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Units.hh"
#include "nuchic/ProcessInfo.hh"

// TODO: Turn this into a factory to reduce the number of includes
#include "nuchic/BeamMapper.hh"
#include "nuchic/HadronicMapper.hh"
#include "nuchic/FinalStateMapper.hh"
#include "nuchic/PhaseSpaceMapper.hh"
#include "nuchic/QuasielasticTestMapper.hh"
#include "plugins/Sherpa/Channels1.hh"
#include "plugins/Sherpa/Channels3.hh"

#include "plugins/Sherpa/SherpaMEs.hh"

#include "yaml-cpp/yaml.h"

nuchic::Channel<nuchic::FourVector> BuildChannelTest(const YAML::Node &node, std::shared_ptr<nuchic::Beam> beam) {
    nuchic::Channel<nuchic::FourVector> channel;
    channel.mapping = std::make_unique<nuchic::QuasielasticTestMapper>(node, beam);
    nuchic::AdaptiveMap2 map(channel.mapping -> NDims(), 2);
    channel.integrator = nuchic::Vegas2(map, {});
    return channel;
}

template<typename T, typename ...Args>
nuchic::Channel<nuchic::FourVector> BuildChannel(size_t nlep, size_t nhad,
                                                 std::shared_ptr<nuchic::Mapper<nuchic::FourVector>> beam,
                                                 std::shared_ptr<nuchic::Mapper<nuchic::FourVector>> hadron,
                                                 Args... args) {
    nuchic::Channel<nuchic::FourVector> channel;
    auto finalStateMapper = std::make_unique<T>(std::forward<Args>(args)...);
    channel.mapping = std::make_unique<nuchic::PSMapper>(nlep, nhad, beam, hadron, std::move(finalStateMapper));
    nuchic::AdaptiveMap2 map(channel.mapping -> NDims(), 2);
    channel.integrator = nuchic::Vegas2(map, nuchic::VegasParams{});
    return channel;
}

template<typename T, typename ...Args>
nuchic::Channel<nuchic::FourVector> BuildChannelSherpa(size_t nlep, size_t nhad,
                                                       std::shared_ptr<nuchic::Mapper<nuchic::FourVector>> beam,
                                                       std::shared_ptr<nuchic::Mapper<nuchic::FourVector>> hadron,
                                                       Args... args) {
    nuchic::Channel<nuchic::FourVector> channel;
    auto sherpaMap = std::make_unique<T>(std::forward<Args>(args)...);
    auto finalStateMapper = std::make_unique<nuchic::SherpaMapper>(nlep+nhad-2, std::move(sherpaMap));
    channel.mapping = std::make_unique<nuchic::PSMapper>(nlep, nhad, beam, hadron, std::move(finalStateMapper));
    nuchic::AdaptiveMap2 map(channel.mapping -> NDims(), 2);
    channel.integrator = nuchic::Vegas2(map, nuchic::VegasParams{});
    return channel;
}

nuchic::EventGen::EventGen(const std::string &configFile, SherpaMEs *const _sherpa) :
  runCascade{false}, outputEvents{false}, sherpa(_sherpa) {
    config = YAML::LoadFile(configFile);

    // Setup random number generator
    auto seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    if(config["Initialize"]["seed"])
        seed = config["Initialize"]["seed"].as<unsigned int>();
    spdlog::trace("Seeding generator with: {}", seed);
    Random::Instance().Seed(seed);

    // Load initial states
    beam = std::make_shared<Beam>(config["Beams"].as<Beam>());
    nucleus = std::make_shared<Nucleus>(config["Nucleus"].as<Nucleus>());

    // Event counts
    total_events = config["EventGen"]["TotalEvents"].as<size_t>();
    nevents = 0; // Initialize to zero
    max_batch = config["EventGen"]["MaxBatch"].as<size_t>();

    // Initialize Cascade parameters
    spdlog::debug("Cascade mode: {}", config["Cascade"]["Run"].as<bool>());
    if(config["Cascade"]["Run"].as<bool>()) {
        cascade = std::make_unique<Cascade>(config["Cascade"].as<Cascade>());
    } else {
        cascade = nullptr;
    }

    // Initialize hard cross-sections
    auto scatteringNode = config["Main"]["Hard Scattering"];
    auto runMode = config["Main"]["Run Mode"].as<nuchic::RunMode>();
    scattering = HardScatteringFactory::Create(scatteringNode["Model"].as<std::string>(),
            scatteringNode, runMode);
    scattering -> SetSherpa(sherpa);
    if(runMode == RunMode::FixedAngle)
        scattering -> SetScatteringAngle(config["Main"]["Angle"].as<double>()*1.0_deg);
    else if(runMode == RunMode::FixedAngleEnergy) {
        scattering -> SetScatteringAngle(config["Main"]["Angle"].as<double>()*1.0_deg);
        scattering -> SetFinalLeptonEnergy(config["Main"]["ELepFinal"].as<double>());
    }

    // Initialize the leptonic process
    auto leptonicProcesses = config["Leptonic Tensor"].as<std::vector<nuchic::Process_Info>>();
    for(const auto &beam_id : beam -> BeamIDs()) {
        std::vector<PID> incoming = {nuchic::PID::dummyHadron(), beam_id};
        for(auto info : leptonicProcesses) {
            info.m_ids.insert(info.m_ids.begin(), incoming.begin(), incoming.end());
            for(const auto id : info.m_ids)
                spdlog::info("{}", int(id));
            if(!sherpa->InitializeProcess(info)) {
                spdlog::error("Cannot initialize hard process");
                exit(1);
            }
            scattering -> AddProcess(info);
        }
    }

    // Setup channels
    if(config["TestingPS"]) {
        Channel<FourVector> channel = BuildChannelTest(config["TestingPS"], beam);
        integrand.AddChannel(std::move(channel));
    } else {
        auto beamMapper = std::make_shared<BeamMapper>(1, beam);
        auto hadronMapper = std::make_shared<QESpectralMapper>(0);
        if(scattering -> Processes()[0].m_ids.size() == 4) {
            // Channel<FourVector> channel = BuildChannel<TwoBodyMapper>(2, 2, beamMapper, hadronMapper,
            //                                                            0, pow(Constant::mN, 2));
            // integrand.AddChannel(std::move(channel));
            Channel<FourVector> channel0 = BuildChannelSherpa<PHASIC::C1_0>(2, 2, beamMapper, hadronMapper,
                                                                            0, pow(Constant::mN, 2));
            Channel<FourVector> channel1 = BuildChannelSherpa<PHASIC::C1_1>(2, 2, beamMapper, hadronMapper,
                                                                            0, pow(Constant::mN, 2));
            Channel<FourVector> channel2 = BuildChannelSherpa<PHASIC::C1_2>(2, 2, beamMapper, hadronMapper,
                                                                            0, pow(Constant::mN, 2));
            integrand.AddChannel(std::move(channel0));
            integrand.AddChannel(std::move(channel1));
            integrand.AddChannel(std::move(channel2));
        } else if(scattering -> Processes()[0].m_ids.size() == 6) {
            spdlog::info("Initializing 2->4");
            constexpr double s5 = (Constant::mN)*(Constant::mN);
            Channel<FourVector> channel0 = BuildChannelSherpa<PHASIC::C3_0>(4, 2, beamMapper, hadronMapper,
                                                                            0, 0, 0, s5);
            Channel<FourVector> channel1 = BuildChannelSherpa<PHASIC::C3_1>(4, 2, beamMapper, hadronMapper,
                                                                            0, 0, 0, s5);
            Channel<FourVector> channel2 = BuildChannelSherpa<PHASIC::C3_2>(4, 2, beamMapper, hadronMapper,
                                                                            0, 0, 0, s5);
            Channel<FourVector> channel3 = BuildChannelSherpa<PHASIC::C3_3>(4, 2, beamMapper, hadronMapper,
                                                                            0, 0, 0, s5);
            Channel<FourVector> channel4 = BuildChannelSherpa<PHASIC::C3_4>(4, 2, beamMapper, hadronMapper,
                                                                            0, 0, 0, s5);
            Channel<FourVector> channel5 = BuildChannelSherpa<PHASIC::C3_5>(4, 2, beamMapper, hadronMapper,
                                                                            0, 0, 0, s5);
            Channel<FourVector> channel6 = BuildChannelSherpa<PHASIC::C3_6>(4, 2, beamMapper, hadronMapper,
                                                                            0, 0, 0, s5);
            Channel<FourVector> channel7 = BuildChannelSherpa<PHASIC::C3_7>(4, 2, beamMapper, hadronMapper,
                                                                            0, 0, 0, s5);
            integrand.AddChannel(std::move(channel0));
            integrand.AddChannel(std::move(channel1));
            integrand.AddChannel(std::move(channel2));
            integrand.AddChannel(std::move(channel3));
            integrand.AddChannel(std::move(channel4));
            integrand.AddChannel(std::move(channel5));
            integrand.AddChannel(std::move(channel6));
            integrand.AddChannel(std::move(channel7));
        } else {
            const std::string error = fmt::format("Leptonic Tensor can only handle 2->2 and 2->4 processes. "
                                                  "Got a 2->{} process", leptonicProcesses[0].m_ids.size()-2);
            throw std::runtime_error(error);
        }
    }

    // Setup Multichannel integrator
    // auto params = config["Integration"]["Params"].as<MultiChannelParams>();
    integrator = MultiChannel(integrand.NDims(), integrand.NChannels(), {1000, 10, 1e8, 5e-2, 1});

    // Decide whether to rotate events to be measured w.r.t. the lepton plane
    if(config["Main"]["DoRotate"])
        doRotate = config["Main"]["DoRotate"].as<bool>();

    // Setup Cuts
    if(config["Main"]["HardCuts"])
        doHardCuts = config["Main"]["HardCuts"].as<bool>();
    spdlog::info("Apply hard cuts? {}", doHardCuts);
    hard_cuts = config["HardCuts"].as<nuchic::CutCollection>();

    // if(config["Main"]["EventCuts"])
    //     doEventCuts = config["Main"]["EventCuts"].as<bool>();
    // spdlog::info("Apply event cuts? {}", doEventCuts);
    // event_cuts = config["EventCuts"].as<nuchic::CutCollection>();

    // Setup outputs
    auto output = config["Main"]["Output"];
    if(output["Format"].as<std::string>() == "Nuchic") {
        bool zipped = true;
        if(output["Zipped"])
            zipped = output["Zipped"].as<bool>();
        writer = std::make_unique<NuchicWriter>(output["Name"].as<std::string>(), zipped);
    }
    writer -> WriteHeader(configFile);

    hist = Histogram(200, 0.0, 400.0, "energy");
    hist2 = Histogram(100, 0.0, 800.0, "momentum");
    hist3 = Histogram(50, -1.0, 1.0, "angle");
    hist4 = Histogram(200, 0.0, 1000.0, "energy");
    hist5 = Histogram(200, 0.0, 1000.0, "momentum");
    hist6 = Histogram(50, -1.0, 1.0, "angle");
}

void nuchic::EventGen::Initialize() {
    spdlog::info("Initializing integrator.");
    auto func = [&](const std::vector<FourVector> &mom, const double &wgt) {
        return GenerateEvent(mom, wgt);
    };
    integrand.Function() = func;
    integrator.Optimize(integrand);
    integrator.Summary();
}

void nuchic::EventGen::GenerateEvents() {
    // integrator.Clear();
    // integrator.Set(config["EventGen"]);
    outputEvents = true;
    runCascade = config["Cascade"]["Run"].as<bool>();
    integrator.Parameters().ncalls = 100000;
    integrator(integrand);
    auto result = integrator.Summary();
    std::cout << "Integral = "
        << fmt::format("{:^8.5e} +/- {:^8.5e} ({:^8.5e} %)",
                       result.results.back().Mean(), result.results.back().Error(),
                       result.results.back().Error() / result.results.back().Mean()*100) << std::endl;
    // spdlog::info("Starting generating of n >= {} total events", total_events);
    // spdlog::info("Using a maximum of {} total Vegas batches.", max_batch);
    // Run integrator in batches until the desired number of events are found
    // size_t batch_count = 1;
    // while ((nevents < total_events) & (batch_count <= max_batch)){
    //     integrator.Clear();  // Reset integrator for each batch
    //     auto func = [&](const std::vector<double> &x, const double &wgt) {
    //         auto niterations = config["EventGen"]["iterations"].as<double>();
    //         return Calculate(x, wgt/niterations, batch_count);
    //     };
    //     spdlog::info("Running vegas batch number {}", batch_count);
    //     integrator(func);
    //     spdlog::info("Total events so far: {}/{}", nevents, total_events);
    //     batch_count += 1;
    // }
    // if (batch_count >= max_batch){
    //     spdlog::info("Stopping after reaching max batch threshold.");
    // }

    hist.Save(config["HistTest1"].as<std::string>());
    hist2.Save(config["HistTest2"].as<std::string>());
    hist3.Save(config["HistTest3"].as<std::string>());
    hist4.Save(config["HistTest4"].as<std::string>());
    hist5.Save(config["HistTest5"].as<std::string>());
    hist6.Save(config["HistTest6"].as<std::string>());
}

double nuchic::EventGen::GenerateEvent(const std::vector<FourVector> &mom, const double &wgt) {
    // Initialize the event, which generates the nuclear configuration
    // and initializes the beam particle for the event
    Event event(nucleus, mom, wgt);

    // Initialize the particle ids for the processes
    const auto pids = scattering -> Processes()[0].m_ids;
    for(auto &me : event.MatrixElements()) {
        me.inital_state.resize(2);
        me.inital_state[1] = pids[1];
        for(size_t idx = 2; idx < pids.size() - 1; ++idx)
            me.final_state.push_back(pids[idx]);
    }

    // Calculate the hard cross sections and select one for initial state
    spdlog::debug("Calculating cross section");

    // Obtain the leptonic tensor
    auto leptonTensor = scattering -> LeptonicTensor(event.Momentum(), 100);
    spdlog::trace("Leptonic Tensor: {}", leptonTensor);

    // Obtain the hadronic tensor
    auto hadronTensors = scattering -> HadronicTensor(event);
    spdlog::trace("Hadronic Tensor(proton): {}", hadronTensors.first);
    spdlog::trace("Hadronic Tensor(neutron): {}", hadronTensors.second);
    scattering -> CrossSection(event);
    constexpr double alpha = 1.0/137;
    std::array<std::complex<double>, 16> hTensor, lTensor;
    auto ke = event.Momentum()[1];
    auto kep = event.Momentum()[2];
    auto q = ke - kep;
    auto rotMat = q.AlignZ();
    q = q.Rotate(rotMat);

#if DEBUG_TENSORS
    auto pp = event.Momentum()[0];
    auto ppp = event.Momentum()[3];
    auto e = pp.E();
    pp.E() = sqrt(pp.P2() + pow(Constant::mN, 2));
    auto q2 = q;
    q2.E() = q.E() - e + Constant::mN - pp.E();
    ke = ke.Rotate(rotMat);
    kep = kep.Rotate(rotMat);
    pp = pp.Rotate(rotMat);
    ppp = ppp.Rotate(rotMat);
    auto prefactor = alpha*4*M_PI/pow(q.M2(), 2);
    auto prefactor2 = alpha*4*M_PI;
    auto ppmn2 = pp*ppp-pow(Constant::mN, 2);
    for(size_t mu = 0; mu < 4; ++mu) {
        for(size_t nu = 0; nu < 4; ++nu) {
            // if(mu == 0 || nu == 0) hadronTensor[4*mu+nu] = 0;
            // if(mu == 3 || nu == 3) hadronTensor[4*mu+nu] = 0;
            // if((mu == 2 && nu == 1) || (mu == 1 && nu == 2)) hadronTensor[4*mu+nu] = 0;
            hTensor[4*mu+nu] = 2*(pp[mu]*ppp[nu] + pp[nu]*ppp[mu])*prefactor2;
            lTensor[4*mu+nu] = 2*(ke[mu]*kep[nu] + ke[nu]*kep[mu])*prefactor;
        }
        hTensor[4*mu+mu] += mu == 0 ? -2*ppmn2*prefactor2 : 2*ppmn2*prefactor2;
        lTensor[4*mu+mu] += mu == 0 ? -2*ke*kep*prefactor : 2*ke*kep*prefactor;
    }
#endif

    std::complex<double> amp_p{}, amp_n{};
    const double factor = alpha;
    for(size_t mu = 0; mu < 4; ++mu) {
        for(size_t nu = 0; nu < 4; ++nu) {
            const size_t idx = 4*mu + nu;
            if(nu == 3) {
                hadronTensors.first[idx] = q.E()/q.P()*hadronTensors.first[4*mu];
                hadronTensors.second[idx] = q.E()/q.P()*hadronTensors.second[4*mu];
            } else if(mu == 3) {
                hadronTensors.first[idx] = q.E()/q.P()*hadronTensors.first[nu];
                hadronTensors.second[idx] = q.E()/q.P()*hadronTensors.second[nu];
            }
        }
    }
    for(size_t mu = 0; mu < 4; ++mu) {
        for(size_t nu = 0; nu < 4; ++nu) {
            const size_t idx = 4*mu + nu;
            if((mu == 0 && nu != 0) || (nu == 0 && mu != 0)) {
                amp_p -= hadronTensors.first[idx]*leptonTensor[idx]*factor;
                amp_n -= hadronTensors.second[idx]*leptonTensor[idx]*factor;
            } else {
                amp_p += hadronTensors.first[idx]*leptonTensor[idx]*factor;
                amp_n += hadronTensors.second[idx]*leptonTensor[idx]*factor;
            }
        }
    }

#ifdef CHECK_WARD_ID
    std::vector<double> ward(8);
    for(size_t mu = 0; mu < 4; ++mu) {
        for(size_t nu = 0; nu < 4; ++nu) {
            const size_t idx = 4*mu + nu;
            if(mu == 0) {
                ward[nu] += (q[mu]*hadronTensors.first[idx]).real();
            } else {
                ward[nu] -= (q[mu]*hadronTensors.first[idx]).real();
            }
            if(nu == 0) {
                ward[4+mu] += (q[nu]*hadronTensors.first[idx]).real();
            } else {
                ward[4+mu] -= (q[nu]*hadronTensors.first[idx]).real();
            }
        }
    }
    for(auto &w : ward) w /= amp_p.real();
    if(amp.real() != 0) spdlog::info("Ward Identities: {}", ward);
#endif

    double flux = 1.0/(2*event.Momentum()[1].E())/(2*sqrt(event.Momentum()[0].P2() + Constant::mN2));
    double xsec_p = amp_p.real()*Constant::HBARC2*2*M_PI*flux;
    double xsec_n = amp_n.real()*Constant::HBARC2*2*M_PI*flux;
    double defaultxsec{};
    static double minRatio = std::numeric_limits<double>::infinity();
    static double maxRatio = 0;
    for(size_t i = 0; i < event.MatrixElements().size(); ++i) {
        if(event.CurrentNucleus() -> Nucleons()[i].ID() == PID::proton()) {
            if(event.MatrixElement(i).weight != 0) {
                defaultxsec = event.MatrixElement(i).weight;
            }
            // break;
            event.MatrixElement(i).weight = xsec_p;
        } else {
            event.MatrixElement(i).weight = xsec_n;
        }
    }
    if(defaultxsec != 0) {
        double ratio = xsec_p/defaultxsec;
        if(ratio < minRatio) minRatio = ratio;
        if(ratio > maxRatio) maxRatio = ratio;
        // spdlog::info("Default xsec = {}", defaultxsec);
        // spdlog::info("Sherpa + Noemi xsec = {}", xsec);
        // spdlog::info("Ratio = {}, Range = [{}, {}]", ratio, minRatio, maxRatio);
    }
    if(!scattering -> InitializeEvent(event)) {
        return 0;
    }

    spdlog::trace("Event Phase Space:");
    size_t idx = 0;
    for(const auto &momentum : event.Momentum()) {
        spdlog::trace("\t{}: {}", ++idx, momentum);
    }

    spdlog::trace("Leptons:");
    idx = 0;
    for(const auto &particle : event.Leptons()) {
        spdlog::trace("\t{}: {}", ++idx, particle);
    }

    spdlog::trace("Hadrons:");
    idx = 0;
    for(const auto &particle : event.Hadrons()) {
        spdlog::trace("\t{}: {}", ++idx, particle);
    }

    // Perform hard cuts
    if(doHardCuts) {
        spdlog::debug("Making hard cuts");
        if(!MakeCuts(event)) {
            // Short-circuit the evaluation
            // We want Vegas to adapt to avoid these points, i.e.,
            // the integrand should be interpreted as zero in this region
            return 0;
        }
    }

    // Run the cascade if needed
    if(runCascade) {
        spdlog::debug("Runnning cascade");
        cascade -> Evolve(&event);
    } else {
        for(auto & nucleon : event.CurrentNucleus()->Nucleons()) {
            if(nucleon.Status() == ParticleStatus::propagating) {
                nucleon.Status() = ParticleStatus::escaped;
            }
        }
    }

    // Write out events
    if(outputEvents) {
        // Rotate cuts into plane of outgoing electron before writing
        if (doRotate)
            Rotate(event);
        // Perform event-level final cuts before writing
        bool outputCurrentEvent = true;
        // if(doEventCuts){
        //     spdlog::debug("Making event cuts");
        //     outputCurrentEvent = MakeEventCuts(event);
        // }

        if(outputCurrentEvent) {
            // Keep a running total of the number of surviving events
            nevents += 1;
            spdlog::debug("Found event: {}/{}", nevents, total_events);
            event.Finalize();
            writer -> Write(event);
            const auto energy = (Constant::mN - event.Momentum()[0].E());
            const auto calls = static_cast<double>(integrator.Parameters().ncalls);
            hist.Fill(energy, event.Weight()/calls);
            const auto momentum = event.Momentum()[0].P();
            hist2.Fill(momentum, event.Weight()/calls);
            const auto cosTheta = event.Momentum()[2].CosTheta();
            hist3.Fill(cosTheta, event.Weight()/calls);
            const auto energy_lepton = event.Momentum()[2].E();
            hist4.Fill(energy_lepton, event.Weight()/calls);
            const auto energy_hadron = event.Momentum()[3].P();
            hist5.Fill(energy_hadron, event.Weight()/calls);
            const auto cosTheta_hadron = event.Momentum()[3].CosTheta();
            hist6.Fill(cosTheta_hadron, event.Weight()/calls);
        }
    }

    // Always return the weight when the event passes the initial hard cut.
    // Even if events do not survive the final event-level cuts, Vegas should
    // still interpret the integrand as nonzero in this region.
    return event.Weight();
}

bool nuchic::EventGen::MakeCuts(Event &event) {
    return hard_cuts.EvaluateCuts(event.Particles());
}

// TODO: Create Analysis level cuts
/*
bool nuchic::EventGen::MakeEventCuts(Event &event) {
    // Run through all particles in the event
    for (const auto& pair : event_cuts) {
        auto pid = pair.first;
        auto cut = pair.second;
        bool pid_passed = false;
        for (const auto& particle : event.Particles()){
            // Restrict to matching final-state particles
            if(particle.IsFinal() && particle.ID() == pid)
                // Keep: at least one particle (of a given PID) survives the cut
                if(cut(particle.Momentum())){
                    pid_passed = true;
                    break;
                }
        }
        // Reject: no particles (of a given PID) satisfy the cut
        if(!pid_passed)
            return false;
    }
    return true;
}*/

void nuchic::EventGen::Rotate(Event &event) {
    // Isolate the azimuthal angle of the outgoing electron
    double phi = 0.0;
    for(const auto & particle : event.Particles()){
        if(particle.ID() == PID::electron() && particle.IsFinal()){
            phi = particle.Momentum().Phi();
        }
    }
    // Rotate the coordiantes of particles so that all azimuthal angles phi are
    // measured with respect to the leptonic plane
    std::array<double, 9> rotation = {
        cos(phi),  sin(phi), 0,
        -sin(phi), cos(phi), 0,
        0,         0,        1};
    event.Rotate(rotation);
}
