#ifndef DOUBLE_EXPONENTIAL_HH
#define DOUBLE_EXPONENTIAL_HH

#include <algorithm>
#include <array>
#include <cmath>

#include "nuchic/Integrators/BaseIntegrator.hh"
#include "nuchic/Utilities.hh"

namespace nuchic {

namespace Integrators {

class DoubleExponential : public BaseIntegrator {
    public:
        DoubleExponential(bool cache=false) : BaseIntegrator(cache) {
            Initialize();
        }
        DoubleExponential(const FunctionD &func, bool cache=false) : BaseIntegrator(func, cache) {
            Initialize();
        }
        DoubleExponential(const FunctionVD &func, bool cache=false) 
            : BaseIntegrator(func, cache) { Initialize(); }

        double Integrate(const double&, const double&, double&, double&);
        std::vector<double> IntegrateVec(const double&, const double&,
                                         double&, double&);

    private:
        void Initialize() {
            if(!init) {
                init = true;
                GenerateTables();
            }
        }
        // Setup abscissa and weights for Double Exponential integration
        static constexpr size_t _phases = 8;
        static constexpr size_t _size = 6*1<<_phases;

        struct DEPoint {
            double abscissa;
            double weights; 
        };

        static double t2(const size_t &k) {
            return exp(static_cast<double>(k)/ipow(2.0,_phases));
        }

        static double u1(const size_t &k) {
            return 0.5*M_PI/2.0*(t2(k)+1.0/t2(k));
        }

        static double t3(const size_t &k) {
            return exp(0.5*M_PI/2.0*(t2(k)-1.0/t2(k)));   
        }

        static double t4(const size_t &k) {
            return 0.5*(t3(k)+1.0/t3(k));   
        }

        static double GetAbcs(const size_t &k) {
            return 1.0/(t3(k)*t4(k));
        }

        static double GetWeight(const size_t &k) {
            return u1(k)/(t4(k)*t4(k));
        }

        static void GenerateTables(){
            for(size_t i = 0; i < _size; ++i)
                table[i] = { GetAbcs(i), GetWeight(i) }; 

            for(size_t i = 0; i < _phases+1; ++i)
                ip[i] = 1 << i;
        }

        static std::array<DEPoint, _size> table;
        static std::array<size_t, _phases+1> ip;
        static bool init;
};

}

}

#endif
