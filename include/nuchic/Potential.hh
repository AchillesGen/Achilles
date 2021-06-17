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

}

#endif
