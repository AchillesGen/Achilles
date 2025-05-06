#include "Achilles/CascadeInteractions/DeltaInteractions.hh"
#include "Achilles/ClebschGordan.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Distributions.hh"
#include "Achilles/Event.hh"
#include "Achilles/Integrators/DoubleExponential.hh"
#include "Achilles/Interpolation.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/Random.hh"
#include "Achilles/Utilities.hh"
#include <stdexcept>

using namespace achilles;

DeltaInteraction::DeltaInteraction() : 
    sigma_max{{{PID::pionp(), PID::proton()}, 700},
        {{-PID::pionp(), PID::proton()}, 70},
        {{-PID::pionp(), PID::neutron()}, 700},
        {{PID::pionp(), PID::neutron()}, 70},
        {{PID::pion0(), PID::proton()}, 135},
        {{PID::pion0(), PID::neutron()}, 135},
    },
    outgoing{{{PID::pionp(), PID::proton()}, {PID::deltapp()}},
        {{-PID::pionp(), PID::proton()}, {PID::delta0()}},
        {{-PID::pionp(), PID::neutron()}, {PID::deltam()}},
        {{PID::pionp(), PID::neutron()}, {PID::deltap()}},
        {{PID::pion0(), PID::proton()}, {PID::deltap()}},
        {{PID::pion0(), PID::neutron()}, {PID::delta0()}},
    }
{
    decay_handler = DecayHandler{"data/decays.yml"};
    InitializeInterpolators();
    // TestInterpolation();
    // throw;
}

DeltaInteraction::DeltaInteraction(const YAML::Node &node) : DeltaInteraction() {
    node["Mode"].as<DeltaInteraction::Mode>();
}

DeltaInteraction::~DeltaInteraction() {}

std::vector<std::pair<PID, PID>> DeltaInteraction::InitialStates() const {
    return {
        {PID::proton(), PID::proton()},   {PID::proton(), PID::neutron()},
        {PID::neutron(), PID::neutron()}, {PID::proton(), PID::deltam()},
        {PID::deltapp(), PID::neutron()}, {PID::deltap(), PID::neutron()},
        {PID::deltap(), PID::proton()},   {PID::delta0(), PID::proton()},
        {PID::delta0(), PID::neutron()},  {PID::deltapp(), PID::proton()},
        {PID::nstarp(), PID::proton()},   {PID::nstarp(), PID::neutron()},
        {PID::nstar0(), PID::proton()},   {PID::nstar0(), PID::neutron()},
        {PID::pionp(), PID::proton()},    {PID::pionp(), PID::neutron()},
        {PID::pion0(), PID::proton()},    {PID::pion0(), PID::neutron()},
        {-PID::pionp(), PID::proton()},   {-PID::pionp(), PID::neutron()},
        {PID::deltam(), PID::neutron()},
    };
}

InteractionResults DeltaInteraction::CrossSection(Event &event, size_t part1, size_t part2) const {
    const auto &particle1 = event.Hadrons()[part1];
    const auto &particle2 = event.Hadrons()[part2];
    // std::pair<PID, PID> ids = {particle1.ID(), particle2.ID()};

    // double mdelta = ParticleInfo(PID::deltap()).Mass();
    // double wdelta = ParticleInfo(PID::deltap()).Width();
    // double sqrts = (particle1.Momentum() + particle2.Momentum()).M();

    // auto boost = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    // auto q = (particle1.Momentum().Boost(boost)).M();

    // FourVector pdelta{mdelta, 0, 0, 0};
    // double qdelta = 1;
    // double width = pow(q / qdelta, 3) * mdelta / sqrts * pow(vfunc(q) / vfunc(qdelta), 2) *
    // wdelta;

    // auto prop = 0.25 * width * width / (pow(mdelta - sqrts, 2) + 0.25 * width * width);
    // auto sigma = sigma_max.at(ids);

    // return {{outgoing.at(ids), sigma * pow(qdelta / q, 2) * prop}};

    InteractionResults results;
    double sqrts = (particle1.Momentum() + particle2.Momentum()).M() / 1_GeV;
    ThreeVector boostCM = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    double p1CM = particle1.Momentum().Boost(-boostCM).P();

    // NOTE: Particle 2 is always a nucleon in the cascade algorithm
    if(particle1.Info().IsNucleon()) {
        // Elastic component
        double elastic_xsec = NNElastic(sqrts, particle1.ID(), particle2.ID());
        spdlog::debug("NN -> NN: sigma = {}", elastic_xsec);
        results.push_back({{particle1.ID(), particle2.ID()}, elastic_xsec / 2});
        results.push_back({{particle2.ID(), particle1.ID()}, elastic_xsec / 2});

        // Delta component
        // TODO: Clean up this logic
        if(particle1.ID() == PID::neutron() && particle2.ID() == PID::neutron()) {
            PID delta_id = PID::delta0();
            double delta_xsec = SigmaNN2NDeltaInterp(sqrts, p1CM, delta_id);
            results.push_back({{PID::neutron(), delta_id}, delta_xsec / 2});
            results.push_back({{delta_id, PID::neutron()}, delta_xsec / 2});
            spdlog::debug("nn -> nDelta0: sigma = {}", delta_xsec);

            delta_id = PID::deltam();
            delta_xsec = SigmaNN2NDeltaInterp(sqrts, p1CM, delta_id);
            results.push_back({{PID::proton(), delta_id}, delta_xsec / 2});
            results.push_back({{delta_id, PID::proton()}, delta_xsec / 2});
            spdlog::debug("nn -> pDelta-: sigma = {}", delta_xsec);
        } else if(particle1.ID() == PID::proton() && particle2.ID() == PID::proton()) {
            PID delta_id = PID::deltap();
            double delta_xsec = SigmaNN2NDeltaInterp(sqrts, p1CM, delta_id);
            results.push_back({{PID::proton(), delta_id}, delta_xsec / 2});
            results.push_back({{delta_id, PID::proton()}, delta_xsec / 2});
            spdlog::debug("pp -> pDelta+: sigma = {}", delta_xsec);

            delta_id = PID::deltapp();
            delta_xsec = SigmaNN2NDeltaInterp(sqrts, p1CM, delta_id);
            results.push_back({{PID::neutron(), delta_id}, delta_xsec / 2});
            results.push_back({{delta_id, PID::neutron()}, delta_xsec / 2});
            spdlog::debug("pp -> nDelta++: sigma = {}", delta_xsec);
        } else {
            PID delta_id = PID::delta0();
            double delta_xsec = SigmaNN2NDeltaInterp(sqrts, p1CM, delta_id);
            results.push_back({{PID::proton(), delta_id}, delta_xsec / 2});
            results.push_back({{delta_id, PID::proton()}, delta_xsec / 2});
            spdlog::debug("np -> pDelta0: sigma = {}", delta_xsec);

            delta_id = PID::deltap();
            delta_xsec = SigmaNN2NDeltaInterp(sqrts, p1CM, delta_id);
            results.push_back({{PID::neutron(), delta_id}, delta_xsec / 2});
            results.push_back({{delta_id, PID::neutron()}, delta_xsec / 2});
            spdlog::debug("np -> nDelta+: sigma = {}", delta_xsec);
        }
    } else if(particle1.Info().IsResonance()) {
        // NDelta -> NN
        double xsec =
            SigmaNDelta2NN(sqrts, p1CM, particle1.ID(), particle2.ID(), particle1.Momentum().M());
        spdlog::debug("NDelta -> NN: sigma = {}", xsec);
        // TODO: Clean up this logic
        if(particle1.ID() == PID::deltapp() && particle2.ID() == PID::neutron()) {
            results.push_back({{PID::proton(), PID::proton()}, xsec});
        } else if(particle1.ID() == PID::deltap() && particle2.ID() == PID::proton()) {
            results.push_back({{PID::proton(), PID::proton()}, xsec});
        } else if(particle1.ID() == PID::deltap() && particle2.ID() == PID::neutron()) {
            results.push_back({{PID::proton(), PID::neutron()}, xsec / 2});
            results.push_back({{PID::neutron(), PID::proton()}, xsec / 2});
        } else if(particle1.ID() == PID::delta0() && particle2.ID() == PID::proton()) {
            results.push_back({{PID::proton(), PID::neutron()}, xsec / 2});
            results.push_back({{PID::neutron(), PID::proton()}, xsec / 2});
        } else if(particle1.ID() == PID::delta0() && particle2.ID() == PID::neutron()) {
            results.push_back({{PID::neutron(), PID::neutron()}, xsec});
        } else if(particle1.ID() == PID::deltam() && particle2.ID() == PID::proton()) {
            results.push_back({{PID::neutron(), PID::neutron()}, xsec});
        }

        // NDelta -> NDelta
        // NOTE: Assume cross section is the same up to isospin
        xsec = SigmaNDelta2NDelta(particle1, particle2, PID::deltap(), PID::neutron());
        for(const auto &delta_id : {PID::deltapp(), PID::deltap(), PID::delta0(), PID::deltam()}) {
            for(const auto &nuc_id : {PID::proton(), PID::neutron()}) {
                double iso = GetIso(
                    particle1.Info().IntCharge() / 3, particle2.Info().IntCharge() / 3,
                    ParticleInfo(nuc_id).IntCharge() / 3, ParticleInfo(delta_id).IntCharge() / 3);
                if(iso == 0) continue;
                spdlog::debug("{}, {} -> {}, {}: sigma = {}", particle1.ID(), particle2.ID(),
                              delta_id, nuc_id, xsec * iso);
                results.push_back({{delta_id, nuc_id}, xsec * iso / 2});
                results.push_back({{nuc_id, delta_id}, xsec * iso / 2});
            }
        }
    } else if(particle1.Info().IsPion()) {
        // Delta channels
        int charge = particle1.Info().IntCharge() + particle2.Info().IntCharge();
        PID delta_id = charge == 6   ? PID::deltapp()
                       : charge == 3 ? PID::deltap()
                       : charge == 0 ? PID::delta0()
                                     : PID::deltam();
        double xsec = SigmaNPi2Delta(particle1, particle2, delta_id);
        results.push_back({{delta_id}, xsec});
        spdlog::debug("NPi -> Delta: sigma = {}", xsec);
    } else {
        throw std::runtime_error("DeltaInteraction: Invalid interaction");
    }

    return results;
}

std::vector<Particle> DeltaInteraction::GenerateMomentum(const Particle &particle1,
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
        ma = GenerateMass(particle1, particle2, out_ids[0], out_ids[1], ran, sqrts);
        spdlog::trace("ma = {}, pid = {}", ma, out_ids[0]);
    } else {
        ma = info_a.Mass();
    }

    ParticleInfo info_b(out_ids[1]);
    double mb;
    if(info_b.IsResonance()) {
        mb = GenerateMass(particle1, particle2, out_ids[1], out_ids[0], ran, sqrts);
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

    return {Particle{out_ids[0], paOut, particle1.Position()},
            Particle{out_ids[1], pbOut, particle2.Position()}};

    // TODO: Implement higher than 2-body interactions
}

std::pair<double, double> DeltaInteraction::TestNNElastic(double sqrts) const {
    return {NNElastic(sqrts / 1_GeV, PID::proton(), PID::proton()),
            NNElastic(sqrts / 1_GeV, PID::neutron(), PID::proton())};
}

double DeltaInteraction::TestDeltaDSigma(bool iresonance, double sqrts, double mdelta) const {
    return DSigmaDM(iresonance, sqrts, mdelta, PID::deltapp());
}

double DeltaInteraction::TestDeltaSigma(double sqrts) const {
    auto sigma = [&](double mass) { return DSigmaDM(false, sqrts, mass, PID::deltapp()); };
    Integrator::DoubleExponential integrator(sigma);
    return integrator.Integrate(0.938 + 0.138, sqrts - 0.938, 1e-6, 1e-4);
}

double DeltaInteraction::TestNPiSigma(const Particle &p1, const Particle &p2, PID res) const {
    return SigmaNPi2Delta(p1, p2, res);
}

double DeltaInteraction::TestNDelta2NDelta(const Particle &p1, const Particle &p2, PID delta_out,
                                           PID nout) const {
    return SigmaNDelta2NDelta(p1, p2, delta_out, nout);
}

double DeltaInteraction::SigmaNPi2Delta(const Particle &p1, const Particle &p2, PID res) const {
    ParticleInfo res_info{res};
    SpinState i1{p1.Info().IsoSpin(), p1.Info().IsoSpinZ()};
    SpinState i2{p2.Info().IsoSpin(), p2.Info().IsoSpinZ()};
    SpinState i3{res_info.IsoSpin(), static_cast<double>(i1.zaxis + i2.zaxis) / 2};

    spdlog::trace("I({}) = {}, I({}) = {}, I({}) = {}", p1.ID(), i1.total / 2, p2.ID(),
                  i2.total / 2, res, i3.total / 2);

    spdlog::trace("Iz({}) = {}, Iz({}) = {}, Iz({}) = {}", p1.ID(), i1.zaxis / 2, p2.ID(),
                  i2.zaxis / 2, res, i3.zaxis / 2);

    double iso_factor = pow(ClebschGordan(i1, i2, i3), 2);
    spdlog::trace("iso = {}", iso_factor);
    if(iso_factor < 1e-8) return 0;

    double spin_factor = static_cast<double>(res_info.NSpins()) /
                         static_cast<double>(p1.Info().NSpins()) /
                         static_cast<double>(p2.Info().NSpins());

    double res_mass = (p1.Momentum() + p2.Momentum()).M() / 1_GeV;
    double pcm2 =
        Pcm2(pow(res_mass, 2), pow(p1.Info().Mass() / 1_GeV, 2), pow(p2.Info().Mass() / 1_GeV, 2));

    spdlog::trace("pcm2 = {}", pcm2);
    if(pcm2 <= 0) return 0;

    double gamma_tot = GetEffectiveWidth(res, res_mass, p1.Mass() / 1_GeV, p2.Mass() / 1_GeV, 1);
    double gamma_in = decay_handler.BranchingRatio(res, {p1.ID(), p2.ID()}) * gamma_tot;

    spdlog::trace("gamma_in = {}, gamma_tot = {}", gamma_in, gamma_tot);
    if(std::abs(gamma_in) < 1e-8 || std::abs(gamma_tot) < 1e-8) return 0;

    // Old Effenberger ansatz from GiBUU
    return iso_factor * spin_factor * 4 * M_PI / pcm2 * pow(res_mass, 2) * gamma_in * gamma_tot /
           (pow(pow(res_mass, 2) - pow(res_info.Mass() / 1_GeV, 2), 2) +
            gamma_tot * gamma_tot * pow(res_info.Mass() / 1_GeV, 2)) *
           Constant::HBARC2 / 1_GeV / 1_GeV;
}

double DeltaInteraction::GetEffectiveWidth(PID id, double mass, double mass1, double mass2,
                                           size_t angular_mom) const {
    double pole_mass = ParticleInfo(id).Mass() / 1_GeV;
    double pole_width = ParticleInfo(id).Width() / 1_GeV;

    double k = sqrt(Pcm2(mass * mass, mass1 * mass1, mass2 * mass2));
    double k0 = sqrt(Pcm2(pole_mass * pole_mass, mass1 * mass1, mass2 * mass2));

    double rho_mass =
        k / mass * pow(distributions::BlattWeisskopf(k / Constant::HBARC * 1_GeV, angular_mom), 2);
    double rho_pole =
        k0 / pole_mass *
        pow(distributions::BlattWeisskopf(k0 / Constant::HBARC * 1_GeV, angular_mom), 2);

    return pole_width * rho_mass / rho_pole;
}

double DeltaInteraction::TestDeltaDSigmaDOmegaDM(double cost, double sqrts, double mdelta,
                                                 PID delta_id) {
    // NOTE: We just assume the maximum possible mass final state.
    //       This will underestimate the cross section, but it should be small enough for now.
    const double mn = ParticleInfo(PID::neutron()).Mass() / 1_GeV;
    const double mpi = ParticleInfo(PID::pionp()).Mass() / 1_GeV;
    if(sqrts < 2 * mn + mpi || sqrts < mdelta + mn) return 0;

    double pin2 = sqrts * sqrts / 4 - mn * mn;
    double pout2 = std::max((sqrts * sqrts - pow(mdelta + mn, 2)) *
                                (sqrts * sqrts - pow(mdelta - mn, 2)) / (4 * sqrts * sqrts),
                            0.0);

    double e1 = sqrt(mn * mn + pin2);
    double e2 = e1;
    double e3 = sqrt(mn * mn + pout2);
    double e4 = sqrt(mdelta * mdelta + pout2);

    double prop = 0;
    if(mdelta > mn + mpi) {
        double mdel = ParticleInfo(delta_id).Mass() / 1_GeV;
        double width = GetEffectiveWidth(delta_id, mdelta, mn, mpi, 1);
        spdlog::trace("width = {}", width);
        prop = 1 / M_PI * mdelta * width /
               (pow(mdelta * mdelta - mdel * mdel, 2) + pow(mdelta * width, 2));
    }

    // Integrate over Omega (does not depend on phi)
    double t = mdelta * mdelta + mn * mn - 2 * e4 * e2 + 2 * cost * sqrt(pin2 * pout2);
    double u = 2 * mn * mn - 2 * e3 * e2 - 2 * cost * sqrt(pin2 * pout2);

    double mat = 1. / (32 * M_PI * sqrts * sqrts) * MatNN2NDelta(t, u, mdelta * mdelta, delta_id) *
                 Constant::HBARC2 / 1_GeV / 1_GeV;
    spdlog::trace("{}, {}, {}", mat, pout2, prop);
    return mat * sqrt(pout2) * 2 * mdelta * prop;
}

// NOTE: All units are in GeV to ensure cross sections have right units
double DeltaInteraction::NNElastic(double sqrts, PID id1, PID id2) const {
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
        } else if(sqrts < 10) {
            // Fit between plab = 6 GeV and sqrts = 10 GeV
            return 481.795 - 298.236 / (sqrts * sqrts) - 172.885 * sqrts + 5.04276 * sqrts * sqrts -
                   0.113744 * pow(sqrts, 3) + 0.0011898 * pow(sqrts, 4) -
                   132.601 * log(sqrts * sqrts) + 69.1034 * pow(log(sqrts * sqrts), 2);
        } else {
            // Use fit from TOTEM (https://arxiv.org/pdf/1712.06153)
            return 11.84 - 1.617 * log(sqrts * sqrts) + 0.1359 * pow(log(sqrts * sqrts), 2);
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
        } else if(sqrts < 10) {
            // Fit between plab = 6 GeV and sqrts = 10 GeV
            return 481.795 - 298.236 / (sqrts * sqrts) - 172.885 * sqrts + 5.04276 * sqrts * sqrts -
                   0.113744 * pow(sqrts, 3) + 0.0011898 * pow(sqrts, 4) -
                   132.601 * log(sqrts * sqrts) + 69.1034 * pow(log(sqrts * sqrts), 2);
        } else {
            // Use fit from TOTEM (https://arxiv.org/pdf/1712.06153)
            return 11.84 - 1.617 * log(sqrts * sqrts) + 0.1359 * pow(log(sqrts * sqrts), 2);
        }
    }
}

double DeltaInteraction::SigmaNN2NDelta(double sqrts, double pcm, PID delta_id) const {
    auto sigma = [&](double mass) { return DSigmaDM(false, sqrts, mass, delta_id); };
    Integrator::DoubleExponential integrator(sigma);

    double isofactor;
    if(delta_id == PID::deltapp() || delta_id == PID::deltam())
        isofactor = 1;
    else
        isofactor = 1.0 / 3.0;

    const double mn =
        (ParticleInfo(PID::proton()).Mass() + ParticleInfo(PID::neutron()).Mass()) / 2 / 1_GeV;
    const double mpi =
        (2 * ParticleInfo(PID::pionp()).Mass() + ParticleInfo(PID::pion0()).Mass()) / 3 / 1_GeV;
    return integrator.Integrate(mn + mpi, sqrts - mn, 1e-6, 1e-4) * isofactor / (pcm / 1_GeV);
}

double DeltaInteraction::SigmaNDelta2NN(double sqrts, double pcm, PID delta_id, PID nucleon,
                                        double mass) const {
    double isofactor = 0;
    if(delta_id == PID::deltapp() && nucleon == PID::neutron())
        isofactor = 3.0 / 8.0;
    else if(delta_id == PID::deltap())
        isofactor = nucleon == PID::neutron() ? 1.0 / 2.0 : 1.0 / 8.0;
    else if(delta_id == PID::delta0())
        isofactor = nucleon == PID::neutron() ? 1.0 / 8.0 : 1.0 / 2.0;
    else if(delta_id == PID::deltam() && nucleon == PID::proton())
        isofactor = 3.0 / 8.0;

    return DSigmaDM(true, sqrts, mass / 1_GeV, delta_id) * 8.0 / 3.0 * isofactor / (pcm / 1_GeV);
}

double DeltaInteraction::GenerateMass(const Particle &p1, const Particle &p2, PID res, PID other,
                                      Random &ran, double sqrts) const {
    // TODO: Expand to other resonances
    if(res == PID::deltapp() || res == PID::deltap() || res == PID::delta0() ||
       res == PID::deltam()) {
        // TODO: Figure out a better way to handle this,
        // we always choose the heavier particles to ensure that it is always kinematically allowed
        const double mn = ParticleInfo(PID::neutron()).Mass();
        const double mpi = ParticleInfo(PID::pionp()).Mass();
        double smin = pow(mn + mpi, 2);
        double smax = pow(sqrts - ParticleInfo(other).Mass(), 2);

        // Parameters for generating according to BW (should be most efficient)
        double m2 = pow(ParticleInfo(res).Mass(), 2);
        double mw = ParticleInfo(res).Mass() * ParticleInfo(res).Width();
        double ymax = atan((smin - m2) / mw);
        double ymin = atan((smax - m2) / mw);

        bool is_nd_nd = p1.Info().IsResonance();
        double max_val = 0;
        double md = sqrt(m2) / 1_GeV;
        if(is_nd_nd) {
            auto func = [&](double mass) {
                double spectral = SpectralDelta(res, mass / 1_GeV);
                auto val = -DSigmaND2ND(sqrts / 1_GeV, p2.Momentum().M() / 1_GeV,
                                        ParticleInfo(other).Mass() / 1_GeV,
                                        p1.Momentum().M() / 1_GeV, mass / 1_GeV, spectral);
                return val;
            };
            Brent brent(func);
            double m_max = brent.Minimize(sqrt(smin), sqrt(smax));
            double spectral = SpectralDelta(res, m_max / 1_GeV);
            max_val = DSigmaND2ND(sqrts / 1_GeV, p2.Momentum().M() / 1_GeV,
                                  ParticleInfo(other).Mass() / 1_GeV, p1.Momentum().M() / 1_GeV,
                                  m_max / 1_GeV, spectral);
        } else {
            max_val = DSigmaDM(0, sqrts / 1_GeV, md, res);
        }

        while(true) {
            double s = m2 + mw * tan(ymin + ran.Uniform(0.0, 1.0) * (ymax - ymin));
            double mass = sqrt(s);
            double func_val = 0;
            if(is_nd_nd) {
                double spectral = SpectralDelta(res, mass / 1_GeV);
                func_val = DSigmaND2ND(sqrts / 1_GeV, p2.Momentum().M() / 1_GeV,
                                       ParticleInfo(other).Mass() / 1_GeV,
                                       p1.Momentum().M() / 1_GeV, mass / 1_GeV, spectral);
            } else {
                func_val = DSigmaDM(0, sqrts / 1_GeV, mass / 1_GeV, res);
            }
            spdlog::debug("func_val = {}, max_val = {}", func_val, max_val);
            if(func_val / max_val > ran.Uniform(0.0, 1.0)) return mass;
        }
    } else {
        throw std::runtime_error("DeltaInteraction: Only delta resonances implemented for now!");
    }
}

void DeltaInteraction::InitializeInterpolators() {
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

double DeltaInteraction::SigmaNN2NDeltaInterp(double sqrts, double pcm, PID delta_id) const {
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

double DeltaInteraction::DSigmaDMInterp(bool iresonance, double sqrts, double mdelta,
                                        PID delta_id) const {
    if(sqrts < sqrts_min || sqrts > sqrts_max) return DSigmaDM(iresonance, sqrts, mdelta, delta_id);
    if(mdelta < mass_min || mdelta > mass_max) return DSigmaDM(iresonance, sqrts, mdelta, delta_id);

    if(iresonance) return dsigma_res_ndelta(sqrts, mdelta);
    return dsigma_dm_ndelta(sqrts, mdelta);
}

void DeltaInteraction::TestInterpolation() const {
    for(size_t i = 0; i < 2000; ++i) {
        double sqrts = Random::Instance().Uniform(sqrts_min, sqrts_max);
        double interp = SigmaNN2NDeltaInterp(sqrts, 750, PID::deltapp());
        double exact = SigmaNN2NDelta(sqrts, 750, PID::deltapp());
        spdlog::info("sqrts = {}: |{} - {}| = {} ({})", sqrts, exact, interp,
                     std::abs(exact - interp), std::abs(exact - interp) / std::abs(exact));
        // for(size_t j = 0; j < 20; ++j) {
        //     double mass = Random::Instance().Uniform(mass_min, mass_max);
        //     double interp_dm = DSigmaDMInterp(false, sqrts, mass, PID::deltapp());
        //     double exact_dm = DSigmaDM(false, sqrts, mass, PID::deltapp());
        //     spdlog::info("sqrts = {}, mass = {}: |{} - {}| = {} ({})", sqrts, mass,
        //                  exact_dm, interp_dm, std::abs(exact_dm-interp_dm),
        //                  std::abs(exact_dm-interp_dm)/std::abs(exact_dm));
        // }
    }
}

double DeltaInteraction::DSigmaDM(bool iresonance, double sqrts, double mdelta,
                                  PID delta_id) const {
    // TODO: Figure out a better way to handle this,
    // we always choose the heavier particles to ensure that it is always kinematically allowed
    const double mn = ParticleInfo(PID::neutron()).Mass();
    const double mpi = ParticleInfo(PID::pionp()).Mass();
    if(sqrts < 2 * mn + mpi || sqrts < mdelta + mn) return 0;

    double pin2 = sqrts * sqrts / 4 - mn * mn;
    double pout2 = std::max((sqrts * sqrts - pow(mdelta + mn, 2)) *
                                (sqrts * sqrts - pow(mdelta - mn, 2)) / (4 * sqrts * sqrts),
                            0.0);

    double e1 = sqrt(mn * mn + pin2);
    double e2 = e1;
    double e3 = sqrt(mn * mn + pout2);
    double e4 = sqrt(mdelta * mdelta + pout2);

    double prop = 0;
    if(mdelta > mn + mpi) {
        double mdel = ParticleInfo(delta_id).Mass() / 1_GeV;
        double width = GetEffectiveWidth(delta_id, mdelta, mn, mpi, 1);
        prop = 1 / M_PI * mdelta * width /
               (pow(mdelta * mdelta - mdel * mdel, 2) + pow(mdelta * width, 2));
    }

    // Integrate over Omega (does not depend on phi)
    auto dsigmadomega = [&](double cost) {
        double t = mdelta * mdelta + mn * mn - 2 * e4 * e2 + 2 * cost * sqrt(pin2 * pout2);
        double u = 2 * mn * mn - 2 * e3 * e2 - 2 * cost * sqrt(pin2 * pout2);

        double mat = 1. / (32 * M_PI * sqrts * sqrts) *
                     MatNN2NDelta(t, u, mdelta * mdelta, delta_id) * Constant::HBARC2 / 1_GeV /
                     1_GeV;
        if(iresonance) return mat * sqrt(pin2) / 4;
        return mat * sqrt(pout2) * 2 * mdelta * prop;
    };

    Integrator::DoubleExponential integrator(dsigmadomega);
    return integrator.Integrate(-1, 1, 1e-12, 1e-8);
}

double DeltaInteraction::Pcm2(double s, double s1, double s2) const {
    return pow(s + s1 - s2, 2) / (4 * s) - s1;
}

// This calculation is based on the GiBUU implementation of the Dmitriev and Sushkov model
double DeltaInteraction::MatNN2NDelta(double t, double u, double sdelta, PID delta_id) const {
    static constexpr double fps = 2.202, fp = 1.008, lambda2 = 0.63 * 0.63, kappa2 = 0.2 * 0.2,
                            mpi = 0.14;
    static constexpr double mpi2 = mpi * mpi;

    const double mn =
        (ParticleInfo(PID::proton()).Mass() + ParticleInfo(PID::neutron()).Mass()) / 2 / 1_GeV;
    const double mn2 = mn * mn;
    const double gp = fp * 2 * mn / mpi;
    const double mdelta2 = pow(ParticleInfo(delta_id).Mass() / 1_GeV, 2);

    double zt = (Pcm2(mdelta2, mn2, t) + kappa2) / (Pcm2(sdelta, mn2, t) + kappa2);
    double zu = (Pcm2(mdelta2, mn2, u) + kappa2) / (Pcm2(sdelta, mn2, u) + kappa2);

    double mdelta_eff = sqrt(sdelta);

    // Direct terms
    double m1 =
        t * (t - pow(mdelta_eff - mn, 2)) * pow(pow(mdelta_eff + mn, 2) - t, 2) / (3 * sdelta);
    double m2 =
        u * (u - pow(mdelta_eff - mn, 2)) * pow(pow(mdelta_eff + mn, 2) - u, 2) / (3 * sdelta);

    // Interference terms
    double m12 =
        1. / (2 * sdelta) *
        ((t * u + (sdelta - mn2) * (t + u) - pow(sdelta, 2) + pow(mn2, 2)) *
             (t * u + mn * (mdelta_eff + mn) * (sdelta - mn2)) -
         1.0 / 3.0 * (t * u - pow(mdelta_eff + mn, 2) * (t + u) + pow(mdelta_eff + mn, 4)) *
             (t * u - mn * (mdelta_eff - mn) * (sdelta - mn2)));

    // Form factors
    double ft = (lambda2 - mpi2) / (lambda2 - t);
    double fu = (lambda2 - mpi2) / (lambda2 - u);
    double pt = 1 / (t - mpi2);
    double pu = 1 / (u - mpi2);

    return pow(fps * gp / mpi, 2) *
           (pow(ft, 4) * zt * pt * pt * m1 + pow(fu, 4) * zu * pu * pu * m2 +
            pow(ft * fu, 2) * sqrt(zt * zu) * pt * pu * m12);
}

double DeltaInteraction::GetIso(int n1, int delta1, int n2, int delta2) const {
    // Rotate isospins to have incoming proton
    if(n1 == 0) {
        n2 = 1 - n2;
        delta2 = 1 - delta2;
        n1 = 1;
        n2 = 1 - delta1;
    }

    // Take from Table B.8 from O. Buss et al. / Physics Reports 512 (2012) 1–124
    if(delta1 == 2 && delta2 == 2 && n2 == 1)
        return 9.0 / 4.0;
    else if(delta1 == 1) {
        if(delta2 == 2 && n2 == 0)
            return 3;
        else if(delta2 == 1 && n2 == 1)
            return 0.25;
    } else if(delta1 == 0) {
        if(delta2 == 0 && n2 == 1)
            return 0.25;
        else if(delta2 == 1 && n2 == 0)
            return 4;
    } else if(delta1 == -1) {
        if(delta2 == -1 && n2 == 1)
            return 9.0 / 4.0;
        else if(delta2 == 0 && n2 == 0)
            return 3;
    }

    return 0;
}

double DeltaInteraction::SpectralDelta(PID delta, double mass) const {
    double m0 = ParticleInfo(delta).Mass() / 1_GeV;
    double gamma = ParticleInfo(delta).Width() / 1_GeV;

    return 1 / M_PI * mass * gamma / (pow(m0 * m0 - mass * mass, 2) + mass * mass * gamma * gamma);
}

double DeltaInteraction::SigmaNDelta2NDelta(const Particle &p1, const Particle &p2, PID delta_out,
                                            PID nucleon_out) const {
    double sqrts = (p1.Momentum() + p2.Momentum()).M() / 1_GeV;
    double mn2 = ParticleInfo(nucleon_out).Mass() / 1_GeV;

    // Integrate over outgoing mass
    auto dsigmadm = [&](double mu2) {
        double spectral = SpectralDelta(delta_out, mu2);
        return DSigmaND2ND(sqrts, p2.Momentum().M() / 1_GeV, mn2, p1.Momentum().M() / 1_GeV, mu2,
                           spectral);
    };

    Integrator::DoubleExponential integrator(dsigmadm);
    const double mn =
        (ParticleInfo(PID::proton()).Mass() + ParticleInfo(PID::neutron()).Mass()) / 2 / 1_GeV;
    const double mpi =
        (2 * ParticleInfo(PID::pionp()).Mass() + ParticleInfo(PID::pion0()).Mass()) / 3 / 1_GeV;
    return integrator.Integrate(mn + mpi, sqrts - mn, 1e-6, 1e-4);
}

double DeltaInteraction::DSigmaND2ND(double sqrts, double mn1, double mn2, double mu1, double mu2,
                                     double spectral) const {
    static constexpr double eps = 1e-8;
    if(sqrts < mn2 + mu2 + eps) return 0;
    double pcm_initial = sqrt(Pcm2(sqrts * sqrts, mn1 * mn1, mu1 * mu1));
    double pcm_final = sqrt(Pcm2(sqrts * sqrts, mn2 * mn2, mu2 * mu2));
    double factor = pcm_final / pcm_initial * 2 * mu2 / pow(8 * M_PI * sqrts, 2) *
                    Constant::HBARC2 / 1_GeV / 1_GeV;

    double e1 = sqrt(mn1 * mn1 + pcm_initial * pcm_initial);
    double e3 = sqrt(mn2 * mn2 + pcm_final * pcm_final);

    // Integrate over Omega (does not depend on phi)
    auto dsigmadomega = [&](double cost) {
        double sint = sqrt(1 - cost * cost);
        FourVector p1(e1, 0, 0, pcm_initial);
        FourVector p3(e3, pcm_final * sint, 0, pcm_final * cost);
        double t = (p1 - p3).M2();
        return MatNDelta2NDelta(t, mu1, mu2);
    };

    Integrator::DoubleExponential integrator(dsigmadomega);
    return integrator.Integrate(-1, 1, 1e-12, 1e-8) * factor * 2 * M_PI * spectral;
}

// This calculation is based on Eq. (B.11) from O. Buss et al. / Physics Reports 512 (2012) 1–124
double DeltaInteraction::MatNDelta2NDelta(double t, double mu1, double mu2) const {
    static const double mpi =
        (ParticleInfo(PID::pion0()).Mass() + 2 * ParticleInfo(PID::pionp()).Mass()) / 2 / 1_GeV;
    static const double mpi2 = mpi * mpi;

    static const double mn =
        (ParticleInfo(PID::proton()).Mass() + ParticleInfo(PID::neutron()).Mass()) / 2 / 1_GeV;
    static const double mn2 = mn * mn;

    static constexpr double ga = 1.267, fpi = 92.4 / 1_GeV, lambda2 = 0.63 * 0.63;
    const double fnnpi = mpi * ga / (2 * fpi);
    const double fddpi = 9.0 / 5.0 * fnnpi;
    const double ft = (lambda2 - mpi2) / (lambda2 - t);

    return pow(fnnpi * fddpi / mpi2, 2) / 8 * pow(ft, 4) / pow(t - mpi2, 2) * 16 *
           pow(mu1 + mu2, 2) * mn2 * t / (9 * pow(mu1 * mu2, 2)) *
           (-mu1 * mu1 + 2 * mu1 * mu2 - mu2 * mu2 + t) *
           (pow(mu1, 4) - 2 * pow(mu1, 3) * mu2 + 12 * mu1 * mu1 * mu2 * mu2 -
            2 * mu1 * pow(mu2, 3) + pow(mu2, 4) - 2 * mu1 * mu1 * t + 2 * mu1 * mu2 * t -
            2 * mu2 * mu2 * t + t * t);
}
