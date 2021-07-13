#ifndef POTENTIAL_HH
#define POTENTIAL_HH

#include <complex>
#include <memory>
#include <string>
#include <utility>

#include "nuchic/Nucleus.hh"
#include "nuchic/References.hh"

namespace nuchic {

struct PotentialVals {
    double rvector, rscalar;
    double ivector, iscalar;
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
        virtual PotentialVals operator()(const double&, const double&) = 0;
};

class WiringaPotential : public Potential {
    public:
        WiringaPotential(std::shared_ptr<Nucleus> nucleus,
                         const double &rho0) 
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

        PotentialVals operator()(const double &plab, const double &radius) override {
            const double rho = m_nucleus -> Rho(radius);
            const double rho_ratio = rho/m_rho0;
            const double alpha = 15.52*rho_ratio + 24.93*pow(rho_ratio, 2);
            const double beta = -116*rho_ratio;
            const double lambda = 3.29 - 0.373*rho_ratio;

            PotentialVals results{};
            results.rvector = alpha + beta/(1+pow(plab/lambda, 2));
            return results;
        }

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
                    data[8*j+i] = pt1[i*8+j];
                    data[8*(j+8)+i] = pt2[i*8+j];
                    if(j < 6) data[8*(j+16)+i] = pt3[i*6+j];
                }
            }
        }

        std::string GetReference() const override { return m_ref.GetReference(); }

        PotentialVals operator()(const double&, const double&) override;

    private:
        std::shared_ptr<Nucleus> m_nucleus;
        Reference m_ref;
        std::array<double, 22*8> data;

        double Data(size_t i, size_t j) const { return data[8*(i-1) + j-1]; }
        double CalcTerm(double, double, double, double, double) const;
        double CalcTermSurf(double, double, double, double, double) const;

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
