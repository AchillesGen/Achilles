#include "nuchic/NucleonPDF.hh"
#include "nuchic/ParticleInfo.hh"
#include "nuchic/Integrators/DoubleExponential.hh"
#include <numeric>

nuchic::NucleonPDF::NucleonPDF(const std::string &name, int iset) 
    : m_name{name}, m_iset{iset}, m_set{name},
      m_pdf{LHAPDF::mkPDF(m_name, m_iset)},
      m_pdfs{LHAPDF::mkPDFs(m_name)} {}

double nuchic::NucleonPDF::xfxQ2(const nuchic::PID &pid, double x, double q2) const {
    return m_pdf -> xfxQ2(static_cast<int>(pid), x, q2);
}

double nuchic::NucleonPDF::xfxQ2(int pid, double x, double q2) const {
    return m_pdf -> xfxQ2(pid, x, q2);
}

double nuchic::NucleonPDF::xfxQ2(size_t iset, const nuchic::PID &pid, double x, double q2) const {
    return m_pdfs[iset] -> xfxQ2(static_cast<int>(pid), x, q2);
}

double nuchic::NucleonPDF::xfxQ2(size_t iset, int pid, double x, double q2) const {
    return m_pdfs[iset] -> xfxQ2(pid, x, q2);
}

LHAPDF::PDFUncertainty nuchic::NucleonPDF::Uncertainty(const std::vector<double> &data) {
    return m_set.uncertainty(data);
}

double nuchic::NucleonPDF::F2(const nuchic::PID &pid, double x, double q2, size_t order) const {
    if(order > 1)
        throw std::runtime_error("Order greater than 1 is not currently implemented");

    std::array<double, 5> charge2{1.0/9.0, 4.0/9.0, 1.0/9.0, 4.0/9.0, 1.0/9.0};
    if(pid == PID::neutron()) {
        charge2[0] = 4.0/9.0;
        charge2[1] = 1.0/9.0;
    }

    double result = 0; 
    nuchic::Integrators::DoubleExponential integrator; 
    double ceps1 = 1e-8;
    double ceps2 = 1e-6;
    int nf = m_pdf -> alphaS().numFlavorsQ2(q2);
    double mur2 = ipow(2, 2);
    double muf2 = ipow(2, 2);
    for(int i = 1; i <= nf; ++i) {
        double pdfsum = m_pdf -> xfxQ2(i, x, muf2*q2) + m_pdf -> xfxQ2(-i, x, muf2*q2);
        if(order > 0) {
            spdlog::debug("Calculating alphas quark contribution");
            auto funcC2Plus = [&](double z) {
                return C2q1xFPlus(i, x, z, muf2*q2) + C2q1xFPlus(-i, x, z, muf2*q2);
            };
            integrator.SetFunction(funcC2Plus);
            pdfsum += m_pdf -> alphasQ2(mur2*q2)/(2*M_PI)*x*(
                      integrator.Integrate(x, 1, ceps1, ceps2) 
                      - C2q1xF(i, x, muf2*q2) - C2q1xF(-i, x, muf2*q2)); 
            if(muf2 != 1.0) {
                auto funcP1qq = [&](double z) {
                    return P1qqPlus(i, x, z, muf2*q2) + P1qqPlus(-i, x, z, muf2*q2);
                };
                integrator.SetFunction(funcP1qq);
                pdfsum += m_pdf -> alphasQ2(mur2*q2)/(2*M_PI)*x*(
                          integrator.Integrate(x, 1, ceps1, ceps2)
                          - P1qq(i, x, muf2*q2) - P1qq(-i, x, muf2*q2)
                          + 3.0/2.0*m_pdf -> xfxQ2(i, x, muf2*q2)
                          + 3.0/2.0*m_pdf -> xfxQ2(-i, x, muf2*q2)) * log(muf2);
            }
        }
        result += charge2[static_cast<size_t>(i-1)]*(pdfsum);
    }

    if(order > 0) {
        spdlog::debug("Calculating gluon contribution");
        double lambdag = 1.0/(nf)*std::accumulate(charge2.begin(), charge2.end(), 0);
        auto funcC2g = [&](double z) {
            return C2g1xF(x, z, muf2*q2);
        };
        integrator.SetFunction(funcC2g);
        result += m_pdf -> alphasQ2(mur2*q2)/(2*M_PI)*x*lambdag*integrator.Integrate(x, 1, ceps1, ceps2);
        if(muf2 != 1.0) {
            auto funcP1qg = [&](double z) {
                return P1qg(x, z, muf2*q2);
            };
            integrator.SetFunction(funcP1qg);
            result += m_pdf -> alphasQ2(mur2*q2)/(2*M_PI)*x*lambdag*(
                      integrator.Integrate(x, 1, ceps1, ceps2)) * log(muf2);
        }
    }

    return result;
}

double nuchic::NucleonPDF::F1(const nuchic::PID &pid, double x, double q2, size_t order) const {
    return F2(pid, x, q2, order)/(2*x);
}

std::vector<double> nuchic::NucleonPDF::F2All(const nuchic::PID &pid, double x, double q2, size_t order) const {
    if(order != 0)
        throw std::runtime_error("Order greater than 1 is not currently implemented");

    std::array<double, 5> charge2{1.0/9.0, 4.0/9.0, 1.0/9.0, 4.0/9.0, 1.0/9.0};
    if(pid == PID::neutron()) {
        charge2[0] = 4.0/9.0;
        charge2[1] = 1.0/9.0;
    }

    std::vector<double> results(m_pdfs.size());
    for(size_t iset = 0; iset < m_pdfs.size(); ++iset) {
        for(int i = 1; i < 6; ++i) {
            results[iset] += charge2[static_cast<size_t>(i-1)]*(m_pdfs[iset] -> xfxQ2(i, x, q2)
                                                                + m_pdfs[iset] -> xfxQ2(-i, x, q2));
        }
    }

    return results;
}

std::vector<double> nuchic::NucleonPDF::F1All(const nuchic::PID &pid, double x, double q2, size_t order) const {
    auto results = F2All(pid, x, q2, order);
    for(auto &result : results) result /= (2*x);
    return results;
}

double nuchic::NucleonPDF::C2q1xFPlus(int pid, double x, double z, double q2) const {
    double fx1 = m_pdf -> xfxQ2(pid, x/z, q2)/(x/z);
    double fx2 = m_pdf -> xfxQ2(pid, x, q2)/x;
    double gz = (1+z*z)*(log((1-z)/z)-3.0/4.0)+0.25*(1-z)*(9+5*z);
    return (fx1/z-fx2)*CF*gz/(1-z);
}

double nuchic::NucleonPDF::C2q1xF(int pid, double x, double q2) const {
    double fx = m_pdf -> xfxQ2(pid, x, q2)/x; 
    double gx = M_PI*M_PI/3 + 7*x/2.0 + x*x - x*(2+x)/2*log((1+x)/(1-x))-(log(1-x)-3)*log(1-x)-2*PolyLog(2, 1-x);
    return fx*gx;
}

double nuchic::NucleonPDF::C2g1xF(double x, double z, double q2) const {
    int nf = m_pdf -> alphaS().numFlavorsQ2(q2);
    double fx = m_pdf -> xfxQ2(0, x/z, q2)/(x/z);
    double gx = 2*nf*TR*((x*x+(1-x)*(1-x))*log((1-x)/x)-1+8*x*(1-x));
    return fx*gx/z;
}

double nuchic::NucleonPDF::P1qqPlus(int pid, double x, double z, double q2) const {
    double y = x/z;
    double fx1 = m_pdf -> xfxQ2(pid, y, q2)/y*(1+y*y);
    double fx2 = m_pdf -> xfxQ2(pid, x, q2)/x*(1+x*x);
    return CF*(fx1/z - fx2)/(1-z);
}

double nuchic::NucleonPDF::P1qq(int pid, double x, double q2) const {
    return m_pdf -> xfxQ2(pid, x, q2)/x*(1+x*x)*CF*log(x);
}

double nuchic::NucleonPDF::P1qg(double x, double z, double q2) const {
    return m_pdf -> xfxQ2(0, x/z, q2)/(x/z)*TR*((1-x*x)+x*x);
}
