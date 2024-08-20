#include "Achilles/CascadeInteractions/DeltaInteractions.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Distributions.hh"
#include "Achilles/Event.hh"
#include "Achilles/Integrators/DoubleExponential.hh"
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
{}

DeltaInteraction::DeltaInteraction(const YAML::Node &node) : DeltaInteraction() {
    node["Mode"].as<DeltaInteraction::Mode>();
}

std::vector<std::pair<PID, PID>> DeltaInteraction::InitialStates() const {
    return {
        {PID::proton(), PID::proton()},   {PID::proton(), PID::neutron()},
        {PID::deltapp(), PID::neutron()}, {PID::deltap(), PID::neutron()},
        {PID::deltap(), PID::proton()},   {PID::delta0(), PID::proton()},
        {PID::delta0(), PID::neutron()},  {PID::deltam(), PID::proton()},
        {PID::nstarp(), PID::proton()},   {PID::nstarp(), PID::neutron()},
        {PID::nstar0(), PID::proton()},   {PID::nstar0(), PID::neutron()},
        {PID::pionp(), PID::proton()},    {PID::pionp(), PID::neutron()},
        {PID::pion0(), PID::proton()},    {PID::pion0(), PID::neutron()},
        {-PID::pionp(), PID::proton()},   {-PID::pionp(), PID::neutron()},
    };
}

InteractionResults DeltaInteraction::CrossSection(Event &event, size_t part1, size_t part2) const {
    const auto &particle1 = event.Hadrons()[part1];
    const auto &particle2 = event.Hadrons()[part2];
    std::pair<PID, PID> ids = {particle1.ID(), particle2.ID()};

    double mdelta = ParticleInfo(PID::deltap()).Mass();
    double wdelta = ParticleInfo(PID::deltap()).Width();
    double sqrts = (particle1.Momentum() + particle2.Momentum()).M();

    auto boost = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    auto q = (particle1.Momentum().Boost(boost)).M();

    FourVector pdelta{mdelta, 0, 0, 0};
    double qdelta = 1;
    double width = pow(q / qdelta, 3) * mdelta / sqrts * pow(vfunc(q) / vfunc(qdelta), 2) * wdelta;

    auto prop = 0.25 * width * width / (pow(mdelta - sqrts, 2) + 0.25 * width * width);
    auto sigma = sigma_max.at(ids);

    return {{outgoing.at(ids), sigma * pow(qdelta / q, 2) * prop}};
}

std::vector<Particle> DeltaInteraction::GenerateMomentum(const Particle &, const Particle &,
                                                         const std::vector<PID> &, Random &) const {
    return {};
}

std::pair<double, double> DeltaInteraction::TestNNElastic(double sqrts) const {
    return {NNElastic(sqrts / 1_GeV, PID::proton(), PID::proton()),
            NNElastic(sqrts / 1_GeV, PID::neutron(), PID::proton())};
}

double DeltaInteraction::TestInelastic1Pi(double sqrts, size_t type) const {
    return NNInelastic1PiBackground(sqrts, inelastic_1pi_params[type]);
}

double DeltaInteraction::TestDeltaDSigma(bool iresonance, double sqrts, double mdelta) const {
    return DSigmaDM(iresonance, sqrts, mdelta, PID::deltapp());
}

double DeltaInteraction::TestDeltaSigma(double sqrts) const {
    return SigmaNN2NDelta(sqrts, PID::deltapp());
}

double DeltaInteraction::GetEffectiveWidth(PID id, double mass, double mass1, double mass2,
                                           size_t angular_mom) const {
    double pole_mass = ParticleInfo(id).Mass() / 1_GeV;
    double pole_width = ParticleInfo(id).Width() / 1_GeV;

    double k = sqrt(pcm2(mass * mass, mass1 * mass1, mass2 * mass2));
    double k0 = sqrt(pcm2(pole_mass * pole_mass, mass1 * mass1, mass2 * mass2));

    double rho_mass =
        k / mass * pow(distributions::BlattWeisskopf(k / Constant::HBARC * 1_GeV, angular_mom), 2);
    double rho_pole =
        k0 / pole_mass *
        pow(distributions::BlattWeisskopf(k0 / Constant::HBARC * 1_GeV, angular_mom), 2);

    return pole_width * rho_mass / rho_pole;
}

double DeltaInteraction::TestDeltaDSigmaDOmegaDM(double cost, double sqrts, double mdelta,
                                                 PID delta_id) {
    const double mn =
        (ParticleInfo(PID::proton()).Mass() + ParticleInfo(PID::neutron()).Mass()) / 2 / 1_GeV;
    const double mpi =
        (2 * ParticleInfo(PID::pionp()).Mass() + ParticleInfo(PID::pion0()).Mass()) / 3 / 1_GeV;
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
        spdlog::info("width = {}", width);
        prop = 1 / M_PI * mdelta * width /
               (pow(mdelta * mdelta - mdel * mdel, 2) + pow(mdelta * width, 2));
    }

    // Integrate over Omega (does not depend on phi)
    double t = mdelta * mdelta + mn * mn - 2 * e4 * e2 + 2 * cost * sqrt(pin2 * pout2);
    double u = 2 * mn * mn - 2 * e3 * e2 - 2 * cost * sqrt(pin2 * pout2);

    double mat = 1. / (32 * M_PI * sqrts * sqrts) * MatNN2NDelta(t, u, mdelta * mdelta, delta_id) *
                 Constant::HBARC2 / 1_GeV / 1_GeV;
    spdlog::info("{}, {}, {}", mat, pout2, prop);
    return mat * sqrt(pout2) * 2 * mdelta * prop;
}

double DeltaInteraction::NNInelastic1PiBackground(double sqrts,
                                                  const NNInelastic1PiParams &params) const {
    double x =
        (sqrts - 2 * ParticleInfo(PID::proton()).Mass() - ParticleInfo(PID::pion0()).Mass()) /
        5_GeV;
    return params.norm * pow(x, params.n1) * exp(-(params.a * pow(x, params.n2) + params.b * x));
}

// NOTE: All units are in GeV to ensure cross sections have right units
double DeltaInteraction::NNElastic(double sqrts, PID id1, PID id2) const {
    bool same_iso = id1 == id2;
    double mn = (ParticleInfo(id1).Mass() + ParticleInfo(id2).Mass()) / 2 / 1.0_GeV;
    double threshold = sqrts * sqrts - 4 * mn * mn;
    double plab = sqrts / (2 * mn) * sqrt(threshold);

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

double DeltaInteraction::NN2NDelta(const Particle &, const Particle &) {
    return 0;
};

double DeltaInteraction::SigmaNN2NDelta(double sqrts, PID delta_id) const {
    auto sigma = [&](double mass) { return DSigmaDM(false, sqrts, mass, delta_id); };
    Integrator::DoubleExponential integrator(sigma);
    return integrator.Integrate(0.938 + 0.138, sqrts - 0.938, 1e-6, 1e-4);
}

double DeltaInteraction::DSigmaDM(bool iresonance, double sqrts, double mdelta,
                                  PID delta_id) const {
    const double mn =
        (ParticleInfo(PID::proton()).Mass() + ParticleInfo(PID::neutron()).Mass()) / 2 / 1_GeV;
    const double mpi =
        (2 * ParticleInfo(PID::pionp()).Mass() + ParticleInfo(PID::pion0()).Mass()) / 3 / 1_GeV;
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
};

double DeltaInteraction::pcm2(double s, double s1, double s2) const {
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

    double zt = (pcm2(mdelta2, mn2, t) + kappa2) / (pcm2(sdelta, mn2, t) + kappa2);
    double zu = (pcm2(mdelta2, mn2, u) + kappa2) / (pcm2(sdelta, mn2, u) + kappa2);

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
