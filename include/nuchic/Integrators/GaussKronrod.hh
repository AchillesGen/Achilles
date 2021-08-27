#ifndef GAUSS_KRONROD_HH
#define GAUSS_KRONROD_HH

#include "nuchic/Integrators/BaseIntegrator.hh"

namespace nuchic {

namespace Integrators {

class GaussKronrod : public BaseIntegrator {
    public:
        GaussKronrod(bool cache=true) : BaseIntegrator(cache) {}
        GaussKronrod(const FunctionD &func, bool cache=true) 
            : BaseIntegrator(func, cache) {}
        GaussKronrod(const FunctionVD &func, bool cache)
            : BaseIntegrator(func, cache) {}

        double Integrate(const double&, const double&, double&, double&);
        std::vector<double> IntegrateVec(const double&, const double&,
                                         double&, double&);

    private:
        // Variables
        static constexpr size_t knots = 7;
        static constexpr std::array<double, knots+1> KronrodWgts = {0.022935322010529,0.063092092629979,
                                                                    0.104790010322250,0.140653259715525,
                                                                    0.169004726639267,0.190350578064785,
                                                                    0.204432940075298,0.209482141084728};
        static constexpr std::array<double, (knots+1)/2> GaussWgts = {0.129484966168870,0.279705391489277,
                                                                      0.381830050505119,0.417959183673469};
        static constexpr std::array<double, knots> absc = {0.991455371120813,0.949107912342759,
                                                           0.864864423359769,0.741531185599394,
                                                           0.586087235467691,0.405845151377397,
                                                           0.207784955007898};

};

}

}

#endif
