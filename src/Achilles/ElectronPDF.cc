#include "Achilles/ElectronPDF.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Utilities.hh"

using achilles::AlphaQED;
using achilles::ElectronPDF;

double AlphaQED::operator()(double t) const {
    double Q2 = t > 0 ? t : -t;
    size_t i = Q2 < 0.3 ? 0 : Q2 < 3.0 ? 1 : Q2 < 100 ? 2 : 3;
    double sig_lep_gg =
        m_alpha0 / (3 * M_PI) *
        (PiGamma(PID::electron(), Q2) + PiGamma(PID::muon(), Q2) + PiGamma(PID::tau(), Q2));
    double sig_had_gg = m_A[i] + m_B[i] * log(1 + m_C[i] * Q2);
    double sig_top_gg = m_alpha0 / (3 * M_PI) * 3 * PiGamma(PID::top(), Q2);
    double sigma_gg = sig_lep_gg + sig_had_gg + sig_top_gg;

    return m_alpha0 / (1 - sigma_gg);
}

double AlphaQED::PiGamma(PID id, double scale) const {
    double mass2 = ipow(ParticleInfo(id).Mass(), 2);
    if(mass2 == 0)
        throw std::runtime_error("AlphaQED: Can't evolve QED coupling with zero fermion masses");
    double mqs = mass2 / scale;
    if(scale == 0) return 0;
    if(4 * mqs < 1e-4)
        return (-5. / 3 - log(mqs));
    else if(4 * mqs <= 1) {
        double beta = sqrt(1 - 4 * mqs);
        return 1. / 3 - (1 + 2 * mqs) * (2 + beta * log((1 - beta) / (1 + beta)));
    } else {
        return 0;
    }
}

ElectronPDF::ElectronPDF(PID bunch, Scheme scheme, int order, double muh, double alpha0)
    : m_scheme{scheme}, m_bunch{bunch}, m_mass{ParticleInfo(bunch).Mass()}, m_scale{muh},
      m_order{order}, m_alpha{alpha0} {
    m_partons.push_back(bunch);

    double L = log(muh * muh / m_mass / m_mass);
    m_exponent = m_alpha(m_mass * m_mass) / M_PI * (L - 1);
    if(m_order > 3)
        throw std::runtime_error("ElectronPDF: Invalid order. Only valid up to order 3");
}

double ElectronPDF::operator()(double x, double Q2) const {
    if(x >= m_xmax || x <= m_xmin) return 0;

    double alpha = m_alpha(Q2);
    double L = 2 * log(sqrt(Q2) / m_mass);
    double beta_e = 2 * alpha / M_PI * (L - 1);
    double eta = 2 * alpha / M_PI * L;

    double beta{}, beta_s{}, beta_h{};
    switch(m_scheme) {
    case Scheme::eta:
        beta = beta_e;
        beta_s = beta_h = eta;
        break;
    case Scheme::mixed:
        beta = beta_s = beta_e;
        beta_h = eta;
        break;
    case Scheme::beta:
        beta = beta_s = beta_h = beta_e;
    }

    double gamma = std::tgamma(1. + beta / 2);

    // Produces collinear bremsstrahling in exponentiated LLA
    double S = exp(-0.5 * Constant::GAMMA_E * beta + 0.375 * beta_s) / gamma * beta / 2;
    double h0 = -0.25 * (1 + x) * beta_h;

    double h1{}, h2{};
    if(m_order >= 2) {
        h1 = -1. / 32 * beta_h * beta_h *
             ((1 + 3 * x * x) / (1 - x) * log(x) + 4 * (1 + x) * log(1 - x) + 5 + x);
    }
    if(m_order >= 3) {
        h2 = -1. / 384 * ipow(beta_h, 3) *
             ((1 + x) * (6 * Li2(x) + 12 * ipow(log(1 - x), 2) - 3 * M_PI * M_PI) +
              1 / (1 - x) *
                  (1.5 * (1 + 8 * x + 3 * x * x) * log(x) + 6 * (x + 5) * (1 - x) * log(1 - x) +
                   12 * (1 + x * x) * log(x) * log(1 - x) - (0.5 + 3.5 * x * x) * log(x) * log(x) +
                   0.25 * (39 - 24 * x - 15 * x * x)));
    }
    double result = x * (S * pow(1 - x, beta / 2 - 1) + (h0 + h1 + h2));
    if(x > m_xlarge) result *= pow(100, beta / 2) / (pow(100, beta / 2) - 1);
    return result;
}
