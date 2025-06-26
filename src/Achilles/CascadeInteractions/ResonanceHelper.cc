#include "Achilles/CascadeInteractions/ResonanceHelper.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Distributions.hh"
#include "Achilles/Integrators/DoubleExponential.hh"
#include "Achilles/Particle.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/Random.hh"
#include "Achilles/Units.hh"
#include "Achilles/Utilities.hh"

double achilles::resonance::BreitWignerSpectral(PID pid, double mass) {
    double m0 = ParticleInfo(pid).Mass() / 1_GeV;
    double gamma = ParticleInfo(pid).Width() / 1_GeV;

    return 1 / M_PI * mass * gamma / (pow(m0 * m0 - mass * mass, 2) + mass * mass * gamma * gamma);
}

double achilles::resonance::GenerateMass(const Particle &p1, const Particle &p2, PID res, PID other,
                                         Random &ran, double sqrts) {
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
                double spectral = BreitWignerSpectral(res, mass / 1_GeV);
                auto val = -DSigmaND2ND(sqrts / 1_GeV, p2.Momentum().M() / 1_GeV,
                                        ParticleInfo(other).Mass() / 1_GeV,
                                        p1.Momentum().M() / 1_GeV, mass / 1_GeV, spectral);
                return val;
            };
            Brent brent(func);
            double m_max = brent.Minimize(sqrt(smin), sqrt(smax));
            double spectral = BreitWignerSpectral(res, m_max / 1_GeV);
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
                double spectral = BreitWignerSpectral(res, mass / 1_GeV);
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

double achilles::resonance::GetEffectiveWidth(PID id, double mass, double mass1, double mass2,
                                              size_t angular_mom) {
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

double achilles::resonance::DSigmaND2ND(double sqrts, double mn1, double mn2, double mu1,
                                        double mu2, double spectral) {
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
double achilles::resonance::MatNDelta2NDelta(double t, double mu1, double mu2) {
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

// This calculation is based on the GiBUU implementation of the Dmitriev and Sushkov model
double achilles::resonance::MatNN2NDelta(double t, double u, double sdelta, PID delta_id) {
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

double achilles::resonance::DSigmaDM(bool iresonance, double sqrts, double mdelta, PID delta_id) {
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
}

double achilles::resonance::Pcm2(double s, double s1, double s2) {
    return pow(s + s1 - s2, 2) / (4 * s) - s1;
}
