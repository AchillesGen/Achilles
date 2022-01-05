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
#include "nuchic/PhaseSpaceBuilder.hh"
#include "nuchic/BeamMapper.hh"
#include "nuchic/HadronicMapper.hh"
#include "nuchic/FinalStateMapper.hh"
#include "nuchic/PhaseSpaceMapper.hh"
#include "nuchic/QuasielasticTestMapper.hh"
#include "plugins/Sherpa/Channels1.hh"
#include "plugins/Sherpa/Channels3.hh"

#include "plugins/Sherpa/SherpaMEs.hh"
#include "plugins/HepMC3/HepMC3EventWriter.hh"

#include "yaml-cpp/yaml.h"

#include "nuchic/ComplexFmt.hh"

nuchic::Channel<nuchic::FourVector> BuildChannelTest(const YAML::Node &node, std::shared_ptr<nuchic::Beam> beam) {
    nuchic::Channel<nuchic::FourVector> channel;
    channel.mapping = std::make_unique<nuchic::QuasielasticTestMapper>(node, beam);
    nuchic::AdaptiveMap2 map(channel.mapping -> NDims(), 2);
    channel.integrator = nuchic::Vegas2(map, {});
    return channel;
}

template<typename T>
nuchic::Channel<nuchic::FourVector> BuildChannel(size_t nlep, size_t nhad,
                                                 std::shared_ptr<nuchic::Beam> beam,
                                                 const std::vector<double> &masses) {
    nuchic::Channel<nuchic::FourVector> channel;
    channel.mapping = nuchic::PSBuilder(nlep, nhad).Beam(beam, 1)
                                                   .Hadron("QESpectral")
                                                   .FinalState(T::Name(), masses).build();
    nuchic::AdaptiveMap2 map(channel.mapping -> NDims(), 2);
    channel.integrator = nuchic::Vegas2(map, nuchic::VegasParams{});
    return channel;
}

template<typename T>
nuchic::Channel<nuchic::FourVector> BuildChannelSherpa(size_t nlep, size_t nhad,
                                                       std::shared_ptr<nuchic::Beam> beam,
                                                       const std::vector<double> &masses) {
    nuchic::Channel<nuchic::FourVector> channel;
    // channel.mapping = nuchic::PSBuilder(nlep, nhad).Beam(beam, 1)
    //                                                .Hadron("QESpectral")
    //                                                .SherpaFinalState(T::Name(), masses).build();
    channel.mapping = nuchic::PSBuilder(nlep, nhad).Beam(beam, 1)
                                                   .Hadron("QESpectral")
                                                   .SherpaFinalState(T::Name(), masses).build();
    nuchic::AdaptiveMap2 map(channel.mapping -> NDims(), 2);
    channel.integrator = nuchic::Vegas2(map, nuchic::VegasParams{});
    return channel;
}

nuchic::EventGen::EventGen(const std::string &configFile, SherpaMEs *const _sherpa) :
  runCascade{false}, outputEvents{false}, sherpa(_sherpa) {
    config = YAML::LoadFile(configFile);

    // Setup random number generator
    auto seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    if(config["Initialize"]["Seed"])
        if(config["Initialize"]["Seed"].as<int>() > 0)
            seed = config["Initialize"]["Seed"].as<unsigned int>();
    spdlog::trace("Seeding generator with: {}", seed);
    Random::Instance().Seed(seed);

    // Load initial states
    beam = std::make_shared<Beam>(config["Beams"].as<Beam>());
    nucleus = std::make_shared<Nucleus>(config["Nucleus"].as<Nucleus>());

    // Initialize Cascade parameters
    spdlog::debug("Cascade mode: {}", config["Cascade"]["Run"].as<bool>());
    if(config["Cascade"]["Run"].as<bool>()) {
        cascade = std::make_unique<Cascade>(config["Cascade"].as<Cascade>());
    } else {
        cascade = nullptr;
    }

    // Initialize hard cross-sections
    spdlog::debug("Initializing hard interaction");
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
    spdlog::debug("Initializing the leptonic current calculation");
    auto leptonicProcesses = config["Leptonic Tensor"].as<std::vector<nuchic::Process_Info>>();
    for(const auto &beam_id : beam -> BeamIDs()) {
        std::vector<PID> incoming = {nuchic::PID::neutron(), beam_id};
        for(auto info : leptonicProcesses) {
            info.m_ids.insert(info.m_ids.begin(), incoming.begin(), incoming.end());
            for(const auto id : info.m_ids)
                spdlog::debug("{}", int(id));
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
            std::vector<double> masses{0, pow(Constant::mN, 2)};
            // std::vector<double> masses{0, pow(Constant::mN*12, 2)};
            // Channel<FourVector> channel = BuildChannel<TwoBodyMapper>(2, 2, beam, masses);
            // integrand.AddChannel(std::move(channel));
            Channel<FourVector> channel0 = BuildChannelSherpa<PHASIC::C1_0>(2, 2, beam, masses);
            Channel<FourVector> channel1 = BuildChannelSherpa<PHASIC::C1_1>(2, 2, beam, masses);
            Channel<FourVector> channel2 = BuildChannelSherpa<PHASIC::C1_2>(2, 2, beam, masses);
            integrand.AddChannel(std::move(channel0));
            integrand.AddChannel(std::move(channel1));
            integrand.AddChannel(std::move(channel2));
        } else if(scattering -> Processes()[0].m_ids.size() == 6) {
            spdlog::info("Initializing 2->4");
            std::vector<double> masses{0, 0, 0, pow(Constant::mN, 2)};
            // std::vector<double> masses{0, 0, 0, pow(Constant::mN*12, 2)};
            Channel<FourVector> channel0 = BuildChannelSherpa<PHASIC::C3_0>(4, 2, beam, masses);
            Channel<FourVector> channel1 = BuildChannelSherpa<PHASIC::C3_1>(4, 2, beam, masses);
            Channel<FourVector> channel2 = BuildChannelSherpa<PHASIC::C3_2>(4, 2, beam, masses);
            Channel<FourVector> channel3 = BuildChannelSherpa<PHASIC::C3_3>(4, 2, beam, masses);
            Channel<FourVector> channel4 = BuildChannelSherpa<PHASIC::C3_4>(4, 2, beam, masses);
            Channel<FourVector> channel5 = BuildChannelSherpa<PHASIC::C3_5>(4, 2, beam, masses);
            Channel<FourVector> channel6 = BuildChannelSherpa<PHASIC::C3_6>(4, 2, beam, masses);
            Channel<FourVector> channel7 = BuildChannelSherpa<PHASIC::C3_7>(4, 2, beam, masses);
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
    integrator = MultiChannel(integrand.NDims(), integrand.NChannels(), {1000, 2});

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
    bool zipped = true;
    if(output["Zipped"])
        zipped = output["Zipped"].as<bool>();
    if(output["Format"].as<std::string>() == "Nuchic") {
        writer = std::make_unique<NuchicWriter>(output["Name"].as<std::string>(), zipped);
    } else if(output["Format"].as<std::string>() == "HepMC3") {
        writer = std::make_unique<HepMC3Writer>(output["Name"].as<std::string>(), zipped);
    }
    writer -> WriteHeader(configFile);

    hist = Histogram(200, 0.0, 400.0, "energy");
    hist2 = Histogram(100, 0.0, 800.0, "momentum");
    hist3 = Histogram(50, -1.0, 1.0, "angle");
    hist4 = Histogram(200, 0.0, 1000.0, "energy");
    hist5 = Histogram(300, 0.0, 300.0, "observable");
    hist6 = Histogram(200, 0.0, 5.0, "wgt");
}

void nuchic::EventGen::Initialize() {
    // TODO: Clean up loading of previous results
    auto func = [&](const std::vector<FourVector> &mom, const double &wgt) {
        return GenerateEvent(mom, wgt);
    };
    try {
        YAML::Node old_results = YAML::LoadFile("results.yml");
        integrator = old_results["Multichannel"].as<MultiChannel>();
        integrand = old_results["Channels"].as<Integrand<FourVector>>();
        YAML::Node results;
        results["Multichannel"] = integrator;
        results["Channels"] = integrand;
        integrand.Function() = func;
    } catch(const YAML::BadFile &e) {
        spdlog::info("Initializing integrator.");
        integrand.Function() = func;
        if(config["Initialize"]["Accuracy"])
            integrator.Parameters().rtol = config["Initialize"]["Accuracy"].as<double>();
        integrator.Optimize(integrand);
        integrator.Summary();

        YAML::Node results;
        results["Multichannel"] = integrator;
        results["Channels"] = integrand;

        std::ofstream fresults("results.yml");
        fresults << results;
        fresults.close();
    }
}

void nuchic::EventGen::GenerateEvents() {
    outputEvents = true;
    runCascade = config["Cascade"]["Run"].as<bool>();
    integrator.Parameters().ncalls = config["Main"]["NEvents"].as<size_t>();
    integrator(integrand);
    fmt::print("\n");
    auto result = integrator.Summary();
    fmt::print("Integral = {:^8.5e} +/- {:^8.5e} ({:^8.5e} %)\n",
               result.results.back().Mean(), result.results.back().Error(),
               result.results.back().Error() / result.results.back().Mean()*100);

    hist.Save(config["HistTest1"].as<std::string>());
    hist2.Save(config["HistTest2"].as<std::string>());
    hist3.Save(config["HistTest3"].as<std::string>());
    hist4.Save(config["HistTest4"].as<std::string>());
    hist5.Save(config["HistTest5"].as<std::string>());
    hist6.Save(config["HistTest6"].as<std::string>());
}

double nuchic::EventGen::GenerateEvent(const std::vector<FourVector> &mom, const double &wgt) {
    if(outputEvents) {
        static size_t ievent = 0;
        constexpr size_t statusUpdate = 10000;
        if(++ievent % statusUpdate == 0) {
            fmt::print("Generated {} / {} events\r", ievent, config["Main"]["NEvents"].as<size_t>());
        }
    }
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
    auto leptonCurrent = scattering -> LeptonicCurrents(event.Momentum(), 100);
    spdlog::trace("Leptonic Current: {}", leptonCurrent);

    // Obtain the hadronic tensor
    static FFInfoMap protonFF, neutronFF, coherentFF;
    if(protonFF.empty() && neutronFF.empty() && coherentFF.empty()) {
        for(const auto &current : leptonCurrent) {
            protonFF[current.first] = scattering -> FormFactors(PID::proton(), current.first);
            neutronFF[current.first] = scattering -> FormFactors(PID::neutron(), current.first);
            coherentFF[current.first] = scattering -> FormFactors(PID::carbon(), current.first);
        }
    }

    auto hadronCurrent = scattering -> HadronicCurrents(event, protonFF, neutronFF, coherentFF);
    spdlog::trace("Hadronic Current(proton): {}", hadronCurrent[0]);
    spdlog::trace("Hadronic Current(neutron): {}", hadronCurrent[1]);
    // spdlog::trace("Hadronic Current(coherent): {}", hadronCurrent[2]);

#ifdef CHECK_TENSOR
    std::array<std::complex<double>, 16> lTensor, lTensorExact, hTensor, hTensorExact;
    auto pp = event.Momentum()[0];
    auto ppp = event.Momentum()[3];
    auto ke = event.Momentum()[1];
    auto kep = event.Momentum()[2];
    auto q = ke - kep;
    auto rotMat = q.AlignZ();
    q = q.Rotate(rotMat);
    pp = pp.Rotate(rotMat);
    ppp = ppp.Rotate(rotMat);
    ke = ke.Rotate(rotMat);
    kep = kep.Rotate(rotMat);

    auto levicivita = [](int i, int j, int k, int l) -> std::complex<double> {
        return (i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l)/12;
    };
    auto metric = [](int i, int j) -> double {
        return i == j ? i == 0 ? 1 : -1 : 0;
    };

    spdlog::info("pp = {}", pp);
    spdlog::info("ppp = {}", ppp);
    spdlog::info("ke = {}", ke);
    spdlog::info("kep = {}", kep);
    spdlog::info("q = {}", q);
    const auto mw = 79824.4;
    const auto gw = 208.5;
    auto prop2 = std::norm(1.0/((q.M2() - mw*mw) + std::complex<double>(0, 1)*mw*gw));
    auto coupling = pow(0.464864*sqrt(2), 2)*prop2;
    auto coupling2 = pow(0.464864*sqrt(2), 2);
    auto lCurrent = leptonCurrent[24];
    auto hCurrent = hadronCurrent[0][24];
    for(size_t i = 0; i < 4; ++i) {
        auto mu = static_cast<int>(i);
        for(size_t j = 0; j < 4; ++j) {
            auto nu = static_cast<int>(j);
            for(size_t m = 0; m < 4; ++m) {
                lTensor[i*4 + j] += lCurrent[m][i]*std::conj(lCurrent[m][j]);
                hTensor[i*4 + j] += hCurrent[m][i]*std::conj(hCurrent[m][j]);
            }
            lTensorExact[i*4 + j] = (ke[i]*kep[j] + ke[j]*kep[i] - metric(mu, nu)*ke*kep)*coupling;
            hTensorExact[i*4 + j] = (pp[i]*ppp[j] + pp[j]*ppp[i] - metric(mu, nu)*pp*ppp)*coupling2;
            for(size_t k = 0; k < 4; ++k) {
                auto alpha = static_cast<int>(k);
                for(size_t l = 0; l < 4; ++l) {
                    auto beta = static_cast<int>(l);
                    double sign1 = 1;
                    double sign2 = 1;
                    if(k > 0) sign1 = -1;
                    if(l > 0) sign2 = -1;
                    lTensorExact[i*4 + j] -= std::complex<double>(0, 1)*levicivita(alpha, mu, beta, nu)*sign1*sign2*ke[k]*kep[l]*coupling;
                    // hTensorExact[i*4 + j] += std::complex<double>(0, 1)*levicivita(alpha, mu, beta, nu)*sign1*sign2*pp[l]*ppp[k]*coupling2;
                }
            }
            spdlog::info("Sherpa[{}, {}] = {:.3e}", i, j, lTensor[i*4 + j]);
            spdlog::info("LExact[{}, {}] = {:.3e}", i, j, lTensorExact[i*4 + j]);
            spdlog::info("Ratio = {:.3e}", lTensor[i*4 +j] / lTensorExact[i*4 + j]);
            // spdlog::info("HExact[{}, {}] = {:.3e}", i, j, hTensorExact[i*4 + j]);
        }
    }

    std::complex<double> amp_n2{};
    for(size_t i = 0; i < 4; ++i) {
        for(size_t j = 0; j < 4; ++j) {
            double sign = 1;
            if((i == 0 && j != 0) || (j == 0 && i != 0)) sign = -1;
            amp_n2 += sign*lTensorExact[i*4+j]*hTensorExact[i*4+j];
        }
    }
#endif
    
    // scattering -> CrossSection(event);
    double amp2_p{}, amp2_n{}; //, amp2_coh{};
    const size_t nspins = 1 << (event.Momentum().size() - 2);
    for(size_t i = 0; i < nspins; ++i) {
        for(size_t j = 0; j < 4; ++j) {
            double sign = 1.0;
            std::complex<double> amp_p{}, amp_n{}; //, amp_coh{};
            for(size_t mu = 0; mu < 4; ++mu) {
                for(const auto &lcurrent : leptonCurrent) {
                    auto boson = lcurrent.first;
                    if(hadronCurrent[0].find(boson) != hadronCurrent[0].end())
                        amp_p += sign*lcurrent.second[i][mu]*hadronCurrent[0][boson][j][mu];
                    if(hadronCurrent[1].find(boson) != hadronCurrent[1].end())
                        amp_n += sign*lcurrent.second[i][mu]*hadronCurrent[1][boson][j][mu];
                    // if(hadronCurrent[2].find(boson) != hadronCurrent[2].end() && j == 0)
                    //     amp_coh += sign*lcurrent.second[i][mu]*hadronCurrent[2][boson][j][mu];
                }
                sign = -1.0;
            }
            amp2_p += std::norm(amp_p);
            amp2_n += std::norm(amp_n);
            // amp2_coh += std::norm(amp_coh);
        }
    }

#ifdef CHECK_WARD_ID
    auto ke = event.Momentum()[1];
    auto kep = event.Momentum()[2];
    auto q = ke - kep;
    auto rotMat = q.AlignZ();
    q = q.Rotate(rotMat);
    for(size_t i = 0; i < 4; ++i) {
        std::vector<std::complex<double>> ward(3);
        double sign = 1;
        for(size_t mu = 0; mu < 4; ++mu) {
            ward[0] += sign*q[mu]*leptonCurrent[22][i][mu];
            ward[1] += sign*q[mu]*hadronCurrent[0][22][i][mu];
            ward[2] += sign*q[mu]*hadronCurrent[1][22][i][mu];
            sign = -1;
        }
        for(auto &w : ward) w /= amp2_p;
        if(amp2_p != 0) {
            spdlog::info("Spin idx: {}", i);
            spdlog::info("Ward Identities (l, p, n): {: .4e}, {: .4e}, {: .4e}",
                         ward[0], ward[1], ward[2]);
        }
    }
#endif
    double spin_avg = 4;
    if(event.MatrixElement(0).inital_state[1].Abs() == PID::nu_electron() ||
       event.MatrixElement(0).inital_state[1].Abs() == PID::nu_muon() ||
       event.MatrixElement(0).inital_state[1].Abs() == PID::nu_tau())
        spin_avg = 2;
    double flux = 1.0/(2*event.Momentum()[1].E())/(2*sqrt(event.Momentum()[0].P2() + Constant::mN2));
    // double flux = 1.0/(4*sqrt(pow(event.Momentum()[0]*event.Momentum()[1], 2) - event.Momentum()[0].M2()*event.Momentum()[1].M2()));
    // TODO: Make this cleaner for coherent scattering
    // double nuclear_spin = 2;
    static constexpr double to_nb = 1e6;
    double xsec_p = amp2_p*Constant::HBARC2/spin_avg*flux*to_nb;
    double xsec_n = amp2_n*Constant::HBARC2/spin_avg*flux*to_nb;
    // double xsec_coh = amp2_coh*Constant::HBARC2/spin_avg*nuclear_spin*flux*to_nb;
    spdlog::debug("proton xsec = {}", xsec_p);
    spdlog::debug("neutron xsec = {}", xsec_n);
    // spdlog::debug("coherent xsec = {}", xsec_coh);
    for(size_t i = 0; i < event.MatrixElements().size(); ++i) {
        event.MatrixElement(i).inital_state[0] = pids[0];
        event.MatrixElement(i).final_state.push_back(pids[pids.size()-1]);
        if(event.CurrentNucleus() -> Nucleons()[i].ID() == PID::proton()) {
            event.MatrixElement(i).weight = xsec_p;
        } else {
            event.MatrixElement(i).weight = xsec_n;
        }
    }
    // event.CoherentXsec() = xsec_coh;
    if(!scattering -> InitializeEvent(event)) {
        if(outputEvents) {
            event.SetMEWeight(0);
            writer -> Write(event);
        }
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
            if(outputEvents) {
                event.SetMEWeight(0);
                writer -> Write(event);
            }
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
                nucleon.Status() = ParticleStatus::final_state;
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
            const auto omega = event.Momentum()[1].E() - event.Momentum()[2].E();
            hist5.Fill(omega, event.Weight()/calls/(2*M_PI));
            hist6.Fill(log10(event.Weight()));
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
