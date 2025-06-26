#include "Achilles/CascadeInteractions/DeltaInteractions.hh"
#include "Achilles/ClebschGordan.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Event.hh"
#include "Achilles/Integrators/DoubleExponential.hh"
#include "Achilles/Interpolation.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/OsetCrossSections.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/Random.hh"
#include "Achilles/ResonanceHelper.hh"
#include "Achilles/Utilities.hh"
#include <stdexcept>

using namespace achilles;

const DeltaInteraction::AbsorptionModes DeltaInteraction::absorption_modes = {
    {{PID::pionp(), PID::proton()}, {{PID::neutron(), {PID::proton(), PID::proton()}}}},
    {{PID::pionp(), PID::neutron()},
     {{PID::neutron(), {PID::proton(), PID::neutron()}},
      {PID::neutron(), {PID::neutron(), PID::proton()}},
      {PID::proton(), {PID::proton(), PID::proton()}}}},
    {{-PID::pionp(), PID::proton()},
     {{PID::neutron(), {PID::neutron(), PID::neutron()}},
      {PID::proton(), {PID::neutron(), PID::proton()}},
      {PID::proton(), {PID::proton(), PID::neutron()}}}},
    {{-PID::pionp(), PID::neutron()}, {{PID::proton(), {PID::neutron(), PID::neutron()}}}},
    {{PID::pion0(), PID::proton()},
     {{PID::proton(), {PID::proton(), PID::proton()}},
      {PID::neutron(), {PID::proton(), PID::neutron()}},
      {PID::neutron(), {PID::neutron(), PID::proton()}}}},
    {{PID::pion0(), PID::neutron()},
     {{PID::neutron(), {PID::neutron(), PID::neutron()}},
      {PID::proton(), {PID::neutron(), PID::proton()}},
      {PID::proton(), {PID::proton(), PID::neutron()}}}},
};

DeltaInteraction::DeltaInteraction()
    : outgoing{
          {{PID::pionp(), PID::proton()}, {PID::deltapp()}},
          {{-PID::pionp(), PID::proton()}, {PID::delta0()}},
          {{-PID::pionp(), PID::neutron()}, {PID::deltam()}},
          {{PID::pionp(), PID::neutron()}, {PID::deltap()}},
          {{PID::pion0(), PID::proton()}, {PID::deltap()}},
          {{PID::pion0(), PID::neutron()}, {PID::delta0()}},
      } {
    decay_handler = DecayHandler{"data/decays.yml"};
}

DeltaInteraction::DeltaInteraction(const YAML::Node &node) : DeltaInteraction() {
    node["Mode"].as<DeltaInteraction::Mode>();
    exp_sup = node["ExpSup"].as<double>();
}

DeltaInteraction::~DeltaInteraction() {}

std::vector<std::pair<PID, PID>> DeltaInteraction::InitialStates() const {
    return {
        {PID::proton(), PID::deltam()},  {PID::deltapp(), PID::neutron()},
        {PID::deltap(), PID::neutron()}, {PID::deltap(), PID::proton()},
        {PID::delta0(), PID::proton()},  {PID::delta0(), PID::neutron()},
        {PID::deltapp(), PID::proton()}, {PID::nstarp(), PID::proton()},
        {PID::nstarp(), PID::neutron()}, {PID::nstar0(), PID::proton()},
        {PID::nstar0(), PID::neutron()}, {PID::pionp(), PID::proton()},
        {PID::pionp(), PID::neutron()},  {PID::pion0(), PID::proton()},
        {PID::pion0(), PID::neutron()},  {-PID::pionp(), PID::proton()},
        {-PID::pionp(), PID::neutron()}, {PID::deltam(), PID::neutron()},
    };
}

// TODO: Separate out the resonance initial states from pion initial states
InteractionResults DeltaInteraction::CrossSection(Event &event, size_t part1, size_t part2) const {
    const auto &particle1 = event.Hadrons()[part1];
    const auto &particle2 = event.Hadrons()[part2];
    InteractionResults results;
    double sqrts = (particle1.Momentum() + particle2.Momentum()).M() / 1_GeV;
    ThreeVector boostCM = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    double p1CM = particle1.Momentum().Boost(-boostCM).P();

    auto rho = event.CurrentNucleus()->ProtonRho(particle2.Position().Magnitude()) +
               event.CurrentNucleus()->NeutronRho(particle2.Position().Magnitude());
    auto suppression = exp(-exp_sup * rho / Constant::rho0);

    spdlog::debug("rho = {}", rho);
    spdlog::debug("exponential suppression = {}", suppression);
    // NOTE: Particle 2 is always a nucleon in the cascade algorithm
    if(particle1.Info().IsResonance()) {
        // NDelta -> NN
        double xsec = suppression * SigmaNDelta2NN(sqrts, p1CM, particle1.ID(), particle2.ID(),
                                                   particle1.Momentum().M());

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

        // s-channel absorption
        if(auto states = AllowedAbsorption(event, part1, part2); states.size() > 0) {
            OsetCrossSection Oset;
            auto abs_xsec = Oset.SChannelAbsCrossSection(event, part1, part2);
            // Nuclear Physics A568 (1994) 855-872 Table 1
            auto opposite_isospin_xsec = (5. / 6.) * abs_xsec;
            auto same_isospin_xsec = (1. / 6.) * abs_xsec;

            auto same_isospin_counter = 0;
            auto opp_isospin_counter = 0;

            absorption_partners.clear();
            for(const auto &state : states) {
                size_t abs_partner_idx = state.first;
                if(abs_partner_idx == SIZE_MAX) continue;
                absorption_partners.push_back(event.Hadrons()[abs_partner_idx]);
                if(event.Hadrons()[abs_partner_idx].ID() == particle2.ID()) {
                    same_isospin_counter += 1;
                } else {
                    opp_isospin_counter += 1;
                }
            }

            for(const auto &state : states) {
                size_t abs_partner_idx = state.first;
                if(abs_partner_idx == SIZE_MAX) continue;

                if(event.Hadrons()[abs_partner_idx].ID() == particle2.ID()) {
                    results.push_back({state.second, same_isospin_xsec / same_isospin_counter});
                } else {
                    results.push_back({state.second, opposite_isospin_xsec / opp_isospin_counter});
                }
                spdlog::debug("Absorption xsec: {} + {} + {} -> {} + {} = {}",
                              event.Hadrons()[part1].ID(), event.Hadrons()[part2].ID(),
                              event.Hadrons()[abs_partner_idx].ID(), state.second[0],
                              state.second[1], opposite_isospin_xsec / opp_isospin_counter);
            }
        }
    } else {
        throw std::runtime_error("DeltaInteraction: Invalid interaction");
    }

    return results;
}

DeltaInteraction::absorption_states DeltaInteraction::AllowedAbsorption(Event &event, size_t part1,
                                                                        size_t part2) const {
    std::vector<std::pair<size_t, std::vector<PID>>> states;
    // Closest is a pair denoting the closest proton and neutron
    auto closest = FindClosest(event, part1, part2);
    auto modes = absorption_modes.at({event.Hadrons()[part1].ID(), event.Hadrons()[part2].ID()});
    for(const auto &abs_mode : modes) {
        if(abs_mode.absorption_partner == PID::proton()) {
            states.push_back({closest.first, abs_mode.outgoing});
        } else {
            states.push_back({closest.second, abs_mode.outgoing});
        }
    }

    return states;
}

std::pair<size_t, size_t> DeltaInteraction::FindClosest(Event &event, size_t part1,
                                                        size_t part2) const {
    // Particle 2 is the incoming nucleon
    const auto &particle2 = event.Hadrons()[part2];

    double shortest_p_distance = 10000; // fm^2
    double shortest_n_distance = 10000; // fm^2
    size_t closest_p_idx = SIZE_MAX, closest_n_idx = SIZE_MAX;

    for(std::size_t i = 0; i < event.Hadrons().size(); ++i) {
        if(i == part1 || i == part2) continue;
        if(!event.Hadrons()[i].IsBackground()) continue;
        auto distance = (event.Hadrons()[i].Position() - particle2.Position()).Magnitude2();
        if(distance < shortest_p_distance && event.Hadrons()[i].ID() == PID::proton()) {
            shortest_p_distance = distance;
            closest_p_idx = i;
        } else if(distance < shortest_n_distance && event.Hadrons()[i].ID() == PID::neutron()) {
            shortest_n_distance = distance;
            closest_n_idx = i;
        }
    }
    return {closest_p_idx, closest_n_idx};
}

// TODO: Implement higher than 2-body interactions
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

    // S-channel absorption via Oset
    if(particle1.Info().IsPion() && out_ids.size() == 2) {
        return HandleAbsorption(particle1, particle2, out_ids, ran);
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

    return {Particle{out_ids[0], paOut, particle1.Position()},
            Particle{out_ids[1], pbOut, particle2.Position()}};
}

std::vector<Particle> DeltaInteraction::HandleAbsorption(const Particle &particle1,
                                                         const Particle &particle2,
                                                         const std::vector<PID> &out_pids,
                                                         Random &random) const {
    Particle absorption_partner;

    auto partner_charge =
        (ParticleInfo(out_pids[0]).IntCharge() + ParticleInfo(out_pids[1]).IntCharge() -
         particle1.Info().IntCharge() - particle2.Info().IntCharge()) /
        3;

    for(const auto &partner : absorption_partners) {
        if(partner.Info().IntCharge() / 3 == partner_charge) {
            absorption_partner = partner;
            break;
        }
    }
    absorption_partner.Status() = ParticleStatus::absorption_partner;

    // Boost to center of mass
    ThreeVector boostCM =
        (particle1.Momentum() + particle2.Momentum() + absorption_partner.Momentum()).BoostVector();
    spdlog::debug("{}: {}", particle1.ID(), particle1.Momentum());
    spdlog::debug("{}: {}", particle2.ID(), particle2.Momentum());
    spdlog::debug("{}: {}", absorption_partner.ID(), absorption_partner.Momentum());

    FourVector p1CM = particle1.Momentum().Boost(-boostCM);
    FourVector p2CM = particle2.Momentum().Boost(-boostCM);
    FourVector p3CM = absorption_partner.Momentum().Boost(-boostCM);

    FourVector pTotalCM = p1CM + p2CM + p3CM;

    auto s = pTotalCM.M2();
    auto sqrts = sqrt(s);

    auto ma = ParticleInfo(out_pids[0]).Mass();
    auto mb = ParticleInfo(out_pids[1]).Mass();

    const double Eacms = sqrts / 2 * (1 + ma * ma / s - mb * mb / s);
    const double Ebcms = sqrts / 2 * (1 + mb * mb / s - ma * ma / s);
    auto lambda = sqrt(pow(s - ma * ma - mb * mb, 2) - 4 * ma * ma * mb * mb);
    auto pfCMS = lambda / 2 / sqrts;

    std::vector<double> rans(2);
    random.Generate(rans);

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

    spdlog::debug("out1: {}, {}", out_pids[0], paOut.Momentum());
    spdlog::debug("out2: {}, {}", out_pids[1], pbOut.Momentum());

    return {Particle{out_pids[0], paOut, particle1.Position()},
            Particle{out_pids[1], pbOut, particle2.Position()}};
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
    double pcm2 = resonance::Pcm2(pow(res_mass, 2), pow(p1.Info().Mass() / 1_GeV, 2),
                                  pow(p2.Info().Mass() / 1_GeV, 2));

    spdlog::trace("pcm2 = {}", pcm2);
    if(pcm2 <= 0) return 0;

    double gamma_tot =
        resonance::GetEffectiveWidth(res, res_mass, p1.Mass() / 1_GeV, p2.Mass() / 1_GeV, 1);
    double gamma_in = decay_handler.BranchingRatio(res, {p1.ID(), p2.ID()}) * gamma_tot;

    spdlog::trace("gamma_in = {}, gamma_tot = {}", gamma_in, gamma_tot);
    if(std::abs(gamma_in) < 1e-8 || std::abs(gamma_tot) < 1e-8) return 0;

    // Old Effenberger ansatz from GiBUU
    return iso_factor * spin_factor * 4 * M_PI / pcm2 * pow(res_mass, 2) * gamma_in * gamma_tot /
           (pow(pow(res_mass, 2) - pow(res_info.Mass() / 1_GeV, 2), 2) +
            gamma_tot * gamma_tot * pow(res_info.Mass() / 1_GeV, 2)) *
           Constant::HBARC2 / 1_GeV / 1_GeV;
}

// Testing functions
double DeltaInteraction::TestDeltaDSigma(bool iresonance, double sqrts, double mdelta) const {
    return resonance::DSigmaDM(iresonance, sqrts, mdelta, PID::deltapp());
}

double DeltaInteraction::TestDeltaSigma(double sqrts) const {
    auto sigma = [&](double mass) {
        return resonance::DSigmaDM(false, sqrts, mass, PID::deltapp());
    };
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
        double width = resonance::GetEffectiveWidth(delta_id, mdelta, mn, mpi, 1);
        spdlog::trace("width = {}", width);
        prop = 1 / M_PI * mdelta * width /
               (pow(mdelta * mdelta - mdel * mdel, 2) + pow(mdelta * width, 2));
    }

    // Integrate over Omega (does not depend on phi)
    double t = mdelta * mdelta + mn * mn - 2 * e4 * e2 + 2 * cost * sqrt(pin2 * pout2);
    double u = 2 * mn * mn - 2 * e3 * e2 - 2 * cost * sqrt(pin2 * pout2);

    double mat = 1. / (32 * M_PI * sqrts * sqrts) *
                 resonance::MatNN2NDelta(t, u, mdelta * mdelta, delta_id) * Constant::HBARC2 /
                 1_GeV / 1_GeV;
    spdlog::trace("{}, {}, {}", mat, pout2, prop);
    return mat * sqrt(pout2) * 2 * mdelta * prop;
}

double DeltaInteraction::DSigmaDMInterp(bool iresonance, double sqrts, double mdelta,
                                        PID delta_id) const {
    if(sqrts < sqrts_min || sqrts > sqrts_max)
        return resonance::DSigmaDM(iresonance, sqrts, mdelta, delta_id);
    if(mdelta < mass_min || mdelta > mass_max)
        return resonance::DSigmaDM(iresonance, sqrts, mdelta, delta_id);

    if(iresonance) return dsigma_res_ndelta(sqrts, mdelta);
    return dsigma_dm_ndelta(sqrts, mdelta);
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

    return resonance::DSigmaDM(true, sqrts, mass / 1_GeV, delta_id) * 8.0 / 3.0 * isofactor /
           (pcm / 1_GeV);
}

double DeltaInteraction::SigmaNDelta2NDelta(const Particle &p1, const Particle &p2, PID delta_out,
                                            PID nucleon_out) const {
    double sqrts = (p1.Momentum() + p2.Momentum()).M() / 1_GeV;
    double mn2 = ParticleInfo(nucleon_out).Mass() / 1_GeV;

    // Integrate over outgoing mass
    auto dsigmadm = [&](double mu2) {
        double spectral = resonance::BreitWignerSpectral(delta_out, mu2);
        return resonance::DSigmaND2ND(sqrts, p2.Momentum().M() / 1_GeV, mn2,
                                      p1.Momentum().M() / 1_GeV, mu2, spectral);
    };

    Integrator::DoubleExponential integrator(dsigmadm);
    // TODO: Figure out a better way to handle this,
    // we always choose the heavier particles to ensure that it is always kinematically allowed
    const double mn = ParticleInfo(PID::neutron()).Mass() / 1_GeV;
    const double mpi = ParticleInfo(PID::pionp()).Mass() / 1_GeV;
    return integrator.Integrate(mn + mpi, sqrts - mn, 1e-6, 1e-4);
}
