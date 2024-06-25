#ifndef ELECTRON_PDF_HH
#define ELECTRON_PDF_HH

#include "Achilles/PDFBase.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/Units.hh"

namespace achilles {

class AlphaQED {
  private:
    static constexpr std::array<double, 4> m_A{0.0, 0.0, 0.00165, 0.00221};
    static constexpr std::array<double, 4> m_B{0.00835, 0.00238, 0.00299, 0.00293};
    static constexpr std::array<double, 4> m_C{1.0, 3.927, 1.0, 1.0};
    double m_alpha0;

    double PiGamma(PID, double) const;

  public:
    AlphaQED(double alpha0) : m_alpha0(std::move(alpha0)) {}

    double operator()(double) const;
    double AqedThomson() const { return m_alpha0; }
};

class ElectronPDF : public PDFBase {
  public:
    enum class Scheme { beta = 0, eta, mixed };
    ElectronPDF(PID, Scheme, int, double, double);
    double operator()(double, double) const override;

  private:
    static constexpr double m_xmin = 1e-6;
    static constexpr double m_xmax = 0.999999;
    static constexpr double m_xlarge = 0.9999;
    static constexpr double m_q2min = 0.25_GeV * 1.0_GeV;
    static constexpr double m_q2max = 1e14_GeV * 1.0_GeV;
    const Scheme m_scheme;
    const PID m_bunch;
    std::vector<PID> m_partons;
    const double m_mass, m_scale;
    const int m_order;
    double m_exponent;
    AlphaQED m_alpha;
};

} // namespace achilles

#endif
