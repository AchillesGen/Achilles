#ifndef POTENTIAL_HH
#define POTENTIAL_HH

#include <complex>
#include <iostream>
#include <memory>
#include <string>
#include <utility>

#include "nuchic/Nucleus.hh"
#include "nuchic/Particle.hh"
#include "nuchic/References.hh"

#ifdef AUTODIFF
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wuninitialized"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wunused-value"
#include "autodiff/forward/dual.hpp"
#include "autodiff/forward/real.hpp"
#pragma GCC diagnostic pop
#endif

namespace nuchic {

template<typename T>
T sech(T x) {
    return 1.0/cosh(x);
}

template<typename T>
struct PotentialVals {
    T rvector{}, rscalar{};
    T ivector{}, iscalar{};
};

class Potential {
    public:
        Potential() = default;
        Potential(const Potential&) = default;
        Potential(Potential&&) = default;
        Potential& operator=(const Potential&) = default;
        Potential& operator=(Potential&&) = default;

        virtual ~Potential() = default;
        virtual std::string GetReference() const = 0;
        virtual std::string Name() const = 0;
#ifdef AUTODIFF
        virtual PotentialVals<autodiff::dual> operator()(const autodiff::dual&, const autodiff::dual&) const = 0;
#else
        virtual PotentialVals<double> operator()(const double&, const double&) const = 0;
#endif
        PotentialVals<double> derivative_p(double p, double r, double h=step) const {
            auto fp = [&](double x){ return this -> operator()(x, r); };
            return stencil5(fp, p, h);
        }

        PotentialVals<double> derivative_r(double p, double r, double h=step) const {
            auto fr = [&](double x){ return this -> operator()(p, x); };
            return stencil5(fr, r, h);
        }

    private:
        PotentialVals<double> stencil5(std::function<nuchic::PotentialVals<double>(double)> f, double x, double h) const {
            auto fp2h = f(x + 2*h);
            auto fph = f(x + h);
            auto fm2h = f(x - 2*h);
            auto fmh = f(x - h);

            nuchic::PotentialVals<double> results{};
            double den = 12*h;
            results.rvector = (-fp2h.rvector + 8*fph.rvector - 8*fmh.rvector + fm2h.rvector)/den;
            results.ivector = (-fp2h.ivector + 8*fph.ivector - 8*fmh.ivector + fm2h.ivector)/den;
            results.rscalar = (-fp2h.rscalar + 8*fph.rscalar - 8*fmh.rscalar + fm2h.rscalar)/den;
            results.iscalar = (-fp2h.iscalar + 8*fph.iscalar - 8*fmh.iscalar + fm2h.iscalar)/den;
            return results;
        }

        static constexpr double step = 0.01;
};

class WiringaPotential : public Potential {
    public:
        WiringaPotential(std::shared_ptr<Nucleus> nucleus,
                         const double &rho0=0.16) 
                         : m_nucleus{std::move(nucleus)}, m_rho0{rho0}, m_ref{"article", "PhysRevC.38.2967"} {
            m_ref.AddField("title", "{Single-particle potential in dense nuclear matter}");
            m_ref.AddField("author", "{Wiringa, R. B.}");
            m_ref.AddField("journal", "{Phys. Rev. C}");
            m_ref.AddField("volume", "{38}");
            m_ref.AddField("issue", "{6}");
            m_ref.AddField("pages", "{2967--2970}");
            m_ref.AddField("numpages", "{0}");
            m_ref.AddField("year", "{1998}");
            m_ref.AddField("month", "{Dec}");
            m_ref.AddField("publisher", "{American Physical Society}");
            m_ref.AddField("doi", "{10.1103/PhysRevC.38.2967}");
            m_ref.AddField("url", "{https://link.aps.org/doi/10.1103/PhysRevC.38.2967");
        }

        std::string GetReference() const override { return m_ref.GetReference(); }
        std::string Name() const override { return "Wiringa"; }
        double Rho0() const { return m_rho0; }

#ifdef AUTODIFF
        PotentialVals<autodiff::dual> operator()(const autodiff::dual &plab, const autodiff::dual &radius) const override {
            const autodiff::dual rho = m_nucleus -> Rho(static_cast<double>(radius));
            const autodiff::dual rho_ratio = rho/m_rho0;
            //const autodiff::dual alpha = 15.52*rho_ratio + 24.93*pow(rho_ratio, 2);
            const autodiff::dual beta = -116*rho_ratio /plab;
            //const autodiff::dual lambda = (3.29 - 0.373*rho_ratio)*nuchic::Constant::HBARC;

            PotentialVals<autodiff::dual> results{};
            //results.rvector = alpha + beta/(1+pow(plab/lambda, 2));
	    results.rvector = -2 * beta * plab;// * pow(lambda,2) / pow(pow(lambda,2) + pow(plab,2),2);
            return results;
        }
#else
        PotentialVals<double> operator()(const double &plab, const double &radius) const override {
            const double rho = m_nucleus -> Rho(radius);
            const double rho_ratio = rho/m_rho0;
            //const double alpha = 15.52*rho_ratio + 24.93*pow(rho_ratio, 2);
            const double beta = -116*rho_ratio /plab;
	    std::cout << "rho" << rho << " "<< radius << std::endl;
            //const double lambda = (3.29 - 0.373*rho_ratio)*nuchic::Constant::HBARC;

            PotentialVals<double> results{};
            //results.rvector = alpha + beta/(1+pow(plab/lambda, 2));
	    results.rvector = -2 * beta * plab;// * pow(lambda,2) / pow(pow(lambda,2) + pow(plab,2),2);
            return results;
        }
#endif

    private:
        std::shared_ptr<Nucleus> m_nucleus;
        double m_rho0;
        Reference m_ref;
};

class CooperPotential : public Potential {
    public:
        CooperPotential(std::shared_ptr<Nucleus> nucleus) 
                : m_nucleus{std::move(nucleus)}, m_ref{"article", "PhysRevC.80.034605"} {
            m_ref.AddField("title", "{Global Dirac optical potential from helium to lead}");
            m_ref.AddField("author", "{Cooper, E. D. and Hama, S. and Clark, B. C.}");
            m_ref.AddField("journal", "{Phys. Rev. C}");
            m_ref.AddField("volume", "{80}");
            m_ref.AddField("issue", "{3}");
            m_ref.AddField("pages", "{034605}");
            m_ref.AddField("numpages", "{5}");
            m_ref.AddField("year", "{2009}");
            m_ref.AddField("month", "{Sep}");
            m_ref.AddField("publisher", "{American Physical Society}");
            m_ref.AddField("doi", "{10.1103/PhysRevC.80.034605}");
            m_ref.AddField("url", "{https://link.aps.org/doi/10.1103/PhysRevC.80.034605}");

            for(size_t i = 0; i < 8; ++i) {
                for(size_t j = 0; j < 8; ++j) {
                    data[8*j+i] = pt1[j*8+i];
                    data[8*(j+8)+i] = pt2[j*8+i];
                    if(j < 6) data[8*(j+16)+i] = pt3[i*6+j];
                }
            }
        }

        std::string GetReference() const override { return m_ref.GetReference(); }
        std::string Name() const override { return "Cooper"; }

#ifdef AUTODIFF
        PotentialVals<autodiff::dual> operator()(const autodiff::dual &plab, const autodiff::dual &radius) const override {
            return evaluate<autodiff::dual>(plab, radius);
        }
#else
        PotentialVals<double> operator()(const double &plab, const double &radius) const override {
            return evaluate<double>(plab, radius);
        }
#endif

        template<typename T>
        PotentialVals<T> evaluate(const T &plab, const T &radius) const {
            const auto tplab = sqrt(plab*plab + pow(nuchic::Constant::mN, 2)) - nuchic::Constant::mN;
            const auto aa = static_cast<double>(m_nucleus -> NNucleons());
            const auto wt = aa * nuchic::Constant::AMU; 
            const auto ee = tplab;
            const auto acb = cbrt(aa);
            const auto caa = aa / (aa + 20);
            const auto y = caa, y2 = y*y, y3 = y*y2, y4 = y2*y2;
            constexpr auto wp = 1.0072545*nuchic::Constant::AMU;
            const auto el = ee+wp;
            const auto wp2 = wp*wp;
            const auto wt2 = wt*wt;
            const auto pcm = sqrt(wt2*(el*el-wp2)/(wp2+wt2+2.0*wt*el));
            const auto epcm = sqrt(wp2+pcm*pcm);
            const auto etcm = sqrt(wt2+pcm*pcm);
            const auto sr = epcm + etcm;
            const auto e = 1000.0/epcm;
            const auto x = e;
            const auto x2 = x*x;
            const auto x3 = x*x2;
            const auto x4 = x2*x2;
            const auto recv = (etcm / sr);
            const auto recs = (wt / sr);
            constexpr double cv1b = 1.0, cv2b = 1.0, cs1b = 1.0, cs2b = 1.0;
            constexpr double av1b = 0.7, av2b = 0.7, as1b = 0.7, as2b = 0.7;
            const auto sumr = -100. * (Data(1, 1)+Data(1, 2)*x+Data(1, 3)*x2+Data(1, 4)*x3+Data(1, 5)*x4
                                                 +Data(1, 6)*y+Data(1, 7)*y2+Data(1, 8)*y3+Data(20,1)*y4);
            const auto rv1 =   cv1b * (Data(2, 1)+Data(2, 2)*x+Data(2, 3)*x2+Data(2, 4)*x3+Data(2, 5)*x4
                                                 +Data(2, 6)*y+Data(2, 7)*y2+Data(2, 8)*y3+Data(21,1)*y4);
            const auto av1 =   av1b * (Data(3, 1)+Data(3, 2)*x+Data(3, 3)*x2+Data(3, 4)*x3+Data(3, 5)*x4
                                                 +Data(3, 6)*y+Data(3, 7)*y2+Data(3, 8)*y3);
            const auto sumi = -15.0 * (Data(4, 1)+Data(4, 2)*x+Data(4, 3)*x2+Data(4, 4)*x3+Data(4, 5)*x4
                                                 +Data(4, 6)*y+Data(4, 7)*y2+Data(4, 8)*y3+Data(20,2)*y4);
            const auto rv2 =   cv2b * (Data(5, 1)+Data(5, 2)*x+Data(5, 3)*x2+Data(5, 4)*x3+Data(5, 5)*x4
                                                 +Data(5, 6)*y+Data(5, 7)*y2+Data(5, 8)*y3+Data(21,2)*y4);
            const auto av2 =   av2b * (Data(6, 1)+Data(6, 2)*x+Data(6, 3)*x2+Data(6, 4)*x3+Data(6, 5)*x4
                                                 +Data(6, 6)*y+Data(6, 7)*y2+Data(6, 8)*y3);
            const auto diffr = 700. * (Data(7, 1)+Data(7, 2)*x+Data(7, 3)*x2+Data(7, 4)*x3+Data(7, 5)*x4
                                                 +Data(7, 6)*y+Data(7, 7)*y2+Data(7, 8)*y3+Data(20,3)*y4);
            const auto rs1 =   cs1b * (Data(8, 1)+Data(8, 2)*x+Data(8, 3)*x2+Data(8, 4)*x3+Data(8, 5)*x4
                                                 +Data(8, 6)*y+Data(8, 7)*y2+Data(8, 8)*y3+Data(21,3)*y4);
            const auto as1 =   as1b * (Data(9, 1)+Data(9, 2)*x+Data(9, 3)*x2+Data(9, 4)*x3+Data(9, 5)*x4
                                                 +Data(9, 6)*y+Data(9, 7)*y2+Data(9, 8)*y3);
            const auto diffi = -150 * (Data(10, 1)+Data(10, 2)*x+Data(10, 3)*x2+Data(10, 4)*x3+Data(10, 5)*x4
                                                  +Data(10, 6)*y+Data(10, 7)*y2+Data(10, 8)*y3+Data(20,4)*y4);
            const auto rs2 =   cs2b * (Data(11, 1)+Data(11, 2)*x+Data(11, 3)*x2+Data(11, 4)*x3+Data(11, 5)*x4
                                                  +Data(11, 6)*y+Data(11, 7)*y2+Data(11, 8)*y3+Data(21,4)*y4);
            const auto as2 =   as2b * (Data(12, 1)+Data(12, 2)*x+Data(12, 3)*x2+Data(12, 4)*x3+Data(12, 5)*x4
                                                  +Data(12, 6)*y+Data(12, 7)*y2+Data(12, 8)*y3);

            const auto vv = 0.5*(sumr+diffr);
            const auto vs = 0.5*(sumr-diffr);
            const auto wv = 0.5*(sumi+diffi);
            const auto ws = 0.5*(sumi-diffi);

            const auto wv2 = -100.0*(Data(13, 1)+Data(13, 2)*x+Data(13, 3)*x2+Data(13, 4)*x3+Data(13, 5)*x4
                                                +Data(13, 6)*y+Data(13, 7)*y2+Data(13, 8)*y3+Data(20, 5)*y4);
            // const auto rv22 =       (Data(14, 1)+Data(14, 2)*x+Data(14, 3)*x2+Data(14, 4)*x3+Data(14, 5)*x4);
            // const auto av22 =   0.7*(Data(15, 1)+Data(15, 2)*x+Data(15, 3)*x2+Data(15, 4)*x3+Data(15, 5)*x4);

            const auto ws2 =  100.0*(Data(16, 1)+Data(16, 2)*x+Data(16, 3)*x2+Data(16, 4)*x3+Data(16, 5)*x4
                                                +Data(16, 6)*y+Data(16, 7)*y2+Data(16, 8)*y3+Data(20, 6)*y4);
            // const auto rs22 =       (Data(17, 1)+Data(17, 2)*x+Data(17, 3)*x2+Data(17, 4)*x3+Data(17, 5)*x4);
            // const auto as22 =   0.7*(Data(18, 1)+Data(18, 2)*x+Data(18, 3)*x2+Data(18, 4)*x3+Data(18, 5)*x4);

            const auto rv22=rv2;
            const auto av22=av2;
            const auto rs22=rs2;
            const auto as22=as2;

            const T rva1 = CalcTerm<T>(recv*vv, rv1, acb, av1, radius);
            const T rva2 = CalcTerm<T>(recv*wv, rv2, acb, av2, radius) + CalcTermSurf<T>(recv*wv2, rv22, acb, av22, radius);
            const T rsa1 = CalcTerm<T>(recs*vs, rs1, acb, as1, radius);
            const T rsa2 = CalcTerm<T>(recs*ws, rs2, acb, as2, radius) + CalcTermSurf<T>(recs*ws2, rs22, acb, as22, radius);

            return {rva1, rsa1, rva2, rsa2};
        }

    private:
        std::shared_ptr<Nucleus> m_nucleus;
        Reference m_ref;
        std::array<double, 22*8> data{};

        double Data(size_t i, size_t j) const { return data[8*(i-1) + j-1]; }
        template<typename T>
        T CalcTerm(T prefact, T real, T acb, T imag, T radius) const {
            const auto s1 = sech<T>(real*acb/imag);
            const auto s2 = sech<T>(radius/imag);
            const auto prod = s1*s2;
            const auto a = s1 - prod;
            const auto b = s2 - prod;
            return prefact*b/(a+b);
        }
        template<typename T>
        T CalcTermSurf(T prefact, T real, T acb, T imag, T radius) const {
            const auto s1 = sech<T>(real*acb/imag);
            const auto s2 = sech<T>(radius/imag);
            const auto prod = s1*s2;
            const auto a = s1 - prod;
            const auto b = s2 - prod;
            return prefact*a*b/(a+b)/(a+b);
        }

        static constexpr std::array<double, 8*8> pt1{
                     2.313828E+01, -1.233380E+02,  2.432987E+02, -2.092616E+02,
                     6.606549E+01,  5.645371E+00, -1.014396E+01,  5.733192E+00,
                    -6.841731E+00,  3.747771E+01, -6.858889E+01,  5.484194E+01,
                    -1.617735E+01,  9.973483E-01, -1.044247E+00,  4.707822E-01,
                     1.688843E+01, -8.404292E+01,  1.550125E+02, -1.251468E+02,
                     3.722361E+01,  2.912163E+00, -3.560279E+00,  1.645079E+00,
                     2.143528E+02, -9.614898E+02,  1.576657E+03, -1.147987E+03,
                     3.136021E+02,  2.226300E+01, -3.363301E+01,  1.704743E+01,
                     7.852790E+00, -4.551995E+01,  9.525029E+01, -8.644521E+01,
                     2.891801E+01,  5.105014E+00, -7.512490E+00,  3.681682E+00,
                    -8.332038E+00,  4.210331E+01, -7.878047E+01,  6.498050E+01,
                    -1.997590E+01,  4.603199E+00, -7.773754E+00,  4.281708E+00,
                     1.771733E+01, -9.514249E+01,  1.826431E+02, -1.507498E+02,
                     4.563918E+01,  5.566540E+00, -1.001441E+01,  5.551678E+00,
                    -6.557412E+00,  3.596972E+01, -6.571700E+01,  5.233473E+01,
                    -1.534723E+01,  1.007355E+00, -9.669271E-01,  4.066760E-01};
        static constexpr std::array<double, 8*8> pt2{
                     2.183389E+01, -1.085667E+02,  2.012175E+02, -1.637375E+02,
                     4.926579E+01,  2.874114E+00, -3.559180E+00,  1.647730E+00,
                    -5.944315E+01,  3.517880E+02, -6.896441E+02,  6.080210E+02,
                    -2.026092E+02, -3.492561E+01,  5.197784E+01, -2.516421E+01,
                     6.414280E-01, -9.685444E+00,  2.992500E+01, -3.446931E+01,
                     1.361126E+01,  5.209321E+00, -7.563083E+00,  3.671130E+00,
                    -7.989428E+00,  4.318908E+01, -8.862678E+01,  7.860778E+01,
                    -2.582041E+01,  7.390220E+00, -1.207679E+01,  6.360783E+00,
                     1.871682E+01, -1.144552E+02,  2.566304E+02, -2.509193E+02,
                     9.003415E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,
                     1.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
                     0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
                     1.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
                     0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
                     3.713738E+01, -1.952369E+02,  3.842865E+02, -3.379853E+02,
                     1.116151E+02,  0.000000E+00,  0.000000E+00,  0.000000E+00};
        static constexpr std::array<double, 8*6> pt3{
                     1.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
                     0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
                     1.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
                     0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
                     0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
                     0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
                     0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
                     0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
                     0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
                     0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
                     0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
                     0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00};
};

}

#endif
