#ifndef DEUTERON_DIS_HH
#define DEUTERON_DIS_HH

#include "nuchic/Constants.hh"
#include "nuchic/NucleonPDF.hh"
#include "nuchic/SpectralFunction.hh"

#include <string>

namespace nuchic {

class DeuteronDIS {
    public:
        DeuteronDIS(const std::string&, const std::string&, int=0);

        double F2(double, double, double, double) const;
        double F1(double, double, double, double) const;
        std::vector<double> F2All(double, double, double, double) const;
        std::vector<double> F1All(double, double, double, double) const;

    private:
        double qvec(double, double) const;
        double qvec2(double, double) const;
        double omega(double, double) const;
        double omega_t(double, double, double) const;
        double q2t(double, double, double) const;
        double phase_space(double, bool=false) const;
        double xn(double, double, double, double) const;

        SpectralFunction m_spectral;
        NucleonPDF m_pdf;
        static constexpr double mass = 2*Constant::mN/1_GeV;
        static constexpr double mn = Constant::mN/1_GeV;
        static constexpr double mn2 = Constant::mN*Constant::mN/1_GeV/1_GeV;
};

}

#endif
