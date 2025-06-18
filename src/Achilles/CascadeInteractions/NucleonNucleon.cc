#include "Achilles/CascadeInteractions/NucleonNucleon.hh"
#include "Achilles/CascadeInteractions/ResonanceHelper.hh"
#include "Achilles/Event.hh"
#include "Achilles/Integrators/DoubleExponential.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Random.hh"

using namespace achilles;

NucleonNucleon::NucleonNucleon() {
    // TODO: Make sure there is only one instance of decay handler throughout the code
    decay_handler = DecayHandler{"data/decays.yml"};
    InitializeInterpolators();
}

NucleonNucleon::NucleonNucleon(const YAML::Node &node) : NucleonNucleon() {
    mode = node["Mode"].as<NucleonNucleon::Mode>();
    resonance_mode = node["ResonanceMode"].as<NucleonNucleon::ResonanceMode>();
}

std::vector<std::pair<PID, PID>> NucleonNucleon::InitialStates() const {
    return {
        {PID::proton(), PID::proton()},
        {PID::proton(), PID::neutron()},
        {PID::neutron(), PID::neutron()},
    };
}

// TODO: Extend to other resonances and make independent of this class
std::vector<std::pair<PID, PID>>
NucleonNucleon::AllowedResonanceStates(const Particle &particle1, const Particle &particle2) const {
    std::vector<std::pair<PID, PID>> allowed_states;
    int charge = particle1.Info().IntCharge() + particle2.Info().IntCharge();
    if(charge == 6) {
        allowed_states.push_back({PID::deltapp(), PID::neutron()});
        allowed_states.push_back({PID::deltap(), PID::proton()});
    } else if(charge == 3) {
        allowed_states.push_back({PID::deltap(), PID::neutron()});
        allowed_states.push_back({PID::delta0(), PID::proton()});
    } else if(charge == 0) {
        allowed_states.push_back({PID::delta0(), PID::neutron()});
        allowed_states.push_back({PID::deltam(), PID::proton()});
    } else {
        throw std::runtime_error("NucleonNucleon: Invalid charge for resonance states");
    }

    return allowed_states;
}

InteractionResults NucleonNucleon::CrossSection(Event &event, size_t part1, size_t part2) const {
    const auto &particle1 = event.Hadrons()[part1];
    const auto &particle2 = event.Hadrons()[part2];
    InteractionResults results;
    double sqrts = (particle1.Momentum() + particle2.Momentum()).M() / 1_GeV;
    ThreeVector boostCM = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    double p1CM = particle1.Momentum().Boost(-boostCM).P();

    // Elastic component
    double elastic_xsec = NNElastic(sqrts, particle1.ID(), particle2.ID());
    spdlog::debug("NN -> NN: sigma = {}", elastic_xsec);
    results.push_back({{particle1.ID(), particle2.ID()}, elastic_xsec / 2});
    results.push_back({{particle2.ID(), particle1.ID()}, elastic_xsec / 2});

    // TODO: Extend to other resonances
    // Delta component
    auto allowed_resonances = AllowedResonanceStates(particle1, particle2);
    for(const auto &[res_id, nucleon_id] : allowed_resonances) {
        double res_xsec = SigmaNN2NDeltaInterp(sqrts, p1CM, res_id);
        spdlog::debug("NN -> NDelta({}): sigma = {}", res_id, res_xsec);
        results.push_back({{nucleon_id, res_id}, res_xsec / 2});
        results.push_back({{res_id, nucleon_id}, res_xsec / 2});
    }

    return results;
}

// TODO: Implement higher than 2-body interactions
std::vector<Particle> NucleonNucleon::GenerateMomentum(const Particle &particle1,
                                                       const Particle &particle2,
                                                       const std::vector<PID> &out_ids,
                                                       Random &ran) const {
    if(out_ids.size() == 1) {
        auto mom = particle1.Momentum() + particle2.Momentum();
        auto position = ran.Uniform(0.0, 1.0) < 0.5 ? particle1.Position() : particle2.Position();
        if(std::isnan(mom.Momentum()[0])) { spdlog::error("Nan momenutm in Npi -> Delta"); }
        return {Particle{out_ids[0], mom, position, ParticleStatus::propagating}};
    }

    // Boost to center of mass
    ThreeVector boostCM = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    FourVector p1CM = particle1.Momentum().Boost(-boostCM);
    FourVector p2CM = particle2.Momentum().Boost(-boostCM);

    FourVector pTotalCM = p1CM + p2CM;
    auto s = pTotalCM.M2();
    auto sqrts = sqrt(s);

    ParticleInfo info_a(out_ids[0]);
    double ma;
    if(info_a.IsResonance()) {
        ma = resonance::GenerateMass(particle1, particle2, out_ids[0], out_ids[1], ran, sqrts);
        spdlog::trace("ma = {}, pid = {}", ma, out_ids[0]);
    } else {
        ma = info_a.Mass();
    }

    ParticleInfo info_b(out_ids[1]);
    double mb;
    if(info_b.IsResonance()) {
        mb = resonance::GenerateMass(particle1, particle2, out_ids[1], out_ids[0], ran, sqrts);
        spdlog::trace("mb = {}, pid = {}", mb, out_ids[1]);
    } else {
        mb = info_b.Mass();
    }

    const double Eacms = sqrts / 2. * (1. + ma * ma / s - mb * mb / s);
    const double Ebcms = sqrts / 2. * (1. + mb * mb / s - ma * ma / s);
    auto lambda = sqrt(pow(s - ma * ma - mb * mb, 2) - 4. * ma * ma * mb * mb);
    auto pfCMS = lambda / 2. / sqrts;

    std::vector<double> rans(2);
    ran.Generate(rans);

    double cos_cms = 2 * rans[0] - 1;
    double sin_cms = sqrt(1. - cos_cms * cos_cms);
    double phi_cms = 2 * M_PI * rans[1];
    double cosphi_cms = cos(phi_cms);
    double sinphi_cms = sin(phi_cms);

    FourVector paOut = FourVector(Eacms, pfCMS * sin_cms * cosphi_cms, pfCMS * sin_cms * sinphi_cms,
                                  pfCMS * cos_cms);
    FourVector pbOut = FourVector(Ebcms, -pfCMS * sin_cms * cosphi_cms,
                                  -pfCMS * sin_cms * sinphi_cms, -pfCMS * cos_cms);
    paOut = paOut.Boost(boostCM);
    pbOut = pbOut.Boost(boostCM);
    if(std::isnan(paOut.Momentum()[0])) {
        spdlog::error("Nan momenutm for paOut");
        spdlog::error("boost = {}", boostCM);
        spdlog::error("paOut = {}", paOut);
        spdlog::error("{}, {}, {}, {}, {}, {}, {}, {}", out_ids[0], out_ids[1], ma, mb, Eacms,
                      Ebcms, lambda, pfCMS);
        spdlog::error("{}, {}", s, pow(ma + mb, 2));
    }
    if(std::isnan(pbOut.Momentum()[0])) {
        spdlog::error("Nan momenutm for pbOut");
        spdlog::error("boost = {}", boostCM);
        spdlog::error("pbOut = {}", pbOut);
    }

    // Decay resonances if needed
    if(resonance_mode == ResonanceMode::Decay) {
        std::vector<Particle> decays_out;
        if(info_a.IsResonance()) {
            auto decays_a = decay_handler.Decay(Particle{out_ids[0], paOut, particle1.Position()});
            decays_out.insert(decays_out.end(), decays_a.begin(), decays_a.end());
            for(auto &decay : decays_a) { decay.Position() = particle2.Position(); }
        } else {
            decays_out.push_back(Particle{out_ids[0], paOut, particle1.Position()});
        }
        if(info_b.IsResonance()) {
            auto decays_b = decay_handler.Decay(Particle{out_ids[1], pbOut, particle2.Position()});
            for(auto &decay : decays_b) { decay.Position() = particle2.Position(); }
            decays_out.insert(decays_out.end(), decays_b.begin(), decays_b.end());
        } else {
            decays_out.push_back(Particle{out_ids[1], pbOut, particle2.Position()});
        }
        return decays_out;
    }

    return {Particle{out_ids[0], paOut, particle1.Position()},
            Particle{out_ids[1], pbOut, particle2.Position()}};
}

// NOTE: All units are in GeV to ensure cross sections have right units
double NucleonNucleon::NNElastic(double sqrts, PID id1, PID id2) const {
    bool same_iso = id1 == id2;
    double mn = (ParticleInfo(id1).Mass() + ParticleInfo(id2).Mass()) / 2 / 1.0_GeV;
    double threshold = sqrts * sqrts - 4 * mn * mn;
    if(threshold < 0) return 0;
    double plab = sqrts / (2 * mn) * sqrt(threshold);
    spdlog::debug("NNElastic: {}, {}, {}, {}, {}, {}", sqrts, id1, id2, mn, threshold, plab);

    if(same_iso) {
        if(plab < 0.425) {
            return 5.12 * mn / threshold + 1.67;
        } else if(plab < 0.8) {
            return 23.5 + 1000 * pow(plab - 0.7, 4);
        } else if(plab < 2) {
            return 1250 / (plab + 50) - 4 * pow(plab - 1.3, 2);
        } else if(plab < 6) {
            return 77 / (plab + 1.5);
        } else {
            throw std::domain_error("DeltaInteraction: NNElastic energy out of valid region");
        }
    } else {
        if(plab < 0.525) {
            return 17.05 * mn / threshold - 6.83;
        } else if(plab < 0.8) {
            return 33 + 196 * pow(std::abs(plab - 0.95), 2.5);
        } else if(plab < 2) {
            return 31 / sqrt(plab);
        } else if(plab < 6) {
            return 77 / (plab + 1.5);
        } else {
            throw std::domain_error("DeltaInteraction: NNElastic energy out of valid region");
        }
    }
}

double NucleonNucleon::SigmaNN2NDelta(double sqrts, double pcm, PID delta_id) const {
    auto sigma = [&](double mass) { return resonance::DSigmaDM(false, sqrts, mass, delta_id); };
    Integrator::DoubleExponential integrator(sigma);

    double isofactor;
    if(delta_id == PID::deltapp() || delta_id == PID::deltam())
        isofactor = 1;
    else
        isofactor = 1.0 / 3.0;

    return integrator.Integrate(0.938 + 0.138, sqrts - 0.938, 1e-6, 1e-4) * isofactor /
           (pcm / 1_GeV);
}

// Testing functions
std::pair<double, double> NucleonNucleon::TestNNElastic(double sqrts) const {
    return {NNElastic(sqrts / 1_GeV, PID::proton(), PID::proton()),
            NNElastic(sqrts / 1_GeV, PID::neutron(), PID::proton())};
}

void NucleonNucleon::InitializeInterpolators() {
    auto sqrts_vec = Linspace(sqrts_min, sqrts_max, nsqrts);
    // auto mass_vec = Linspace(mass_min, mass_max, nmass);
    std::vector<double> sqrts_valid;
    std::vector<double> dsigma; //, dsigma_dm(nsqrts * nmass), res(nsqrts * nmass);
    for(const auto &sqrts : sqrts_vec) {
        double result = SigmaNN2NDelta(sqrts, 1_GeV, PID::deltapp());
        if(result == 0) continue;
        sqrts_valid.push_back(sqrts);
        dsigma.push_back(result);
        // for(const auto &mass : mass_vec) {
        //     dsigma_dm[idx2] = DSigmaDM(false, sqrts, mass, PID::deltapp());
        //     res[idx2++] = DSigmaDM(true, sqrts, mass, PID::deltapp());
        // }
    }
    sqrts_min = sqrts_valid.front();
    sqrts_max = sqrts_valid.back();

    Interp1D interp_sigma(sqrts_valid, dsigma, InterpolationType::Polynomial);
    interp_sigma.SetPolyOrder(3);
    dsigma_ndelta = interp_sigma;

    // Interp2D interp_dm(sqrts_vec, mass_vec, dsigma_dm);
    // interp_dm.BicubicSpline();
    // dsigma_dm_ndelta = interp_dm;

    // Interp2D interp_res(sqrts_vec, mass_vec, res);
    // interp_res.BicubicSpline();
    // dsigma_res_ndelta = interp_res;
}

double NucleonNucleon::SigmaNN2NDeltaInterp(double sqrts, double pcm, PID delta_id) const {
    if(sqrts < sqrts_min)
        return 0;
    else if(sqrts > sqrts_max) {
        double result = SigmaNN2NDelta(sqrts, pcm, delta_id);
        spdlog::info("Interpolating: {}, {}, {}, {}, {}", sqrts, sqrts_min, sqrts_max, delta_id,
                     result);
        return result;
    }

    double sigma = dsigma_ndelta(sqrts);
    if(sigma < 0) return 0;

    double isofactor;
    if(delta_id == PID::deltapp() || delta_id == PID::deltam())
        isofactor = 1;
    else
        isofactor = 1.0 / 3.0;
    return sigma * isofactor / (pcm / 1_GeV);
}
