#ifndef POTENTIAL_HH
#define POTENTIAL_HH

#include <complex>
#include <iostream>
#include <memory>
#include <string>
#include <utility>

#include "Achilles/Constants.hh"
#include "Achilles/Factory.hh"
#include "Achilles/Particle.hh"
#include "Achilles/References.hh"
#include "Achilles/Utilities.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace achilles {

inline double sech(double x) {
    return 1.0 / cosh(x);
}

struct PotentialVals {
    double rvector{}, rscalar{};
    double ivector{}, iscalar{};
};

class Nucleus;

inline double binomial(size_t n, size_t k) {
    if(k > n)
        throw std::runtime_error(fmt::format("Binomial coefficient requires n > k: {}, {}", n, k));

    double result = 1.0;
    if(k > n - k) k = n - k;
    for(size_t i = 0; i < k; ++i) {
        result *= static_cast<double>(n - i);
        result /= static_cast<double>(i + 1);
    }

    return result;
}

class Potential {
  public:
    Potential() = default;
    Potential(const Potential &) = default;
    Potential(Potential &&) = default;
    Potential &operator=(const Potential &) = default;
    Potential &operator=(Potential &&) = default;

    static std::string Name() { return "Potential"; }

    virtual ~Potential() = default;
    virtual std::string GetReference() const = 0;
    // TODO: Figure out if it should be Nucleus or not passed in
    virtual PotentialVals operator()(double, double) const = 0;

    PotentialVals derivative_p(double p, double r, double h = step) const {
        auto fp = [&](double x) { return this->operator()(x, r); };
        return stencil5(fp, p, h);
    }

    PotentialVals derivative2_p(double p, double r, double h = step) const {
        auto fp = [&](double x) { return this->operator()(x, r); };
        return stencil5second(fp, p, h);
    }

    PotentialVals derivativen_p(size_t n, double p, double r, double h = step) const {
        auto fp = [&](double x) { return this->operator()(x, r); };
        PotentialVals deriv{};
        for(size_t i = 0; i <= n; ++i) {
            auto vals = fp(p + static_cast<double>(i) * h);
            auto sign = pow(-1, static_cast<int>(i + n));
            auto cnk = binomial(n, i);
            deriv.rscalar += sign * cnk * vals.rscalar;
            deriv.iscalar += sign * cnk * vals.iscalar;
            deriv.rvector += sign * cnk * vals.rvector;
            deriv.ivector += sign * cnk * vals.iscalar;
        }

        deriv.rscalar /= pow(h, static_cast<int>(n));
        deriv.iscalar /= pow(h, static_cast<int>(n));
        deriv.rvector /= pow(h, static_cast<int>(n));
        deriv.ivector /= pow(h, static_cast<int>(n));
        return deriv;
    }

    PotentialVals derivative_r(double p, double r, double h = step) const {
        auto fr = [&](double x) { return this->operator()(p, x); };
        return stencil5(fr, r, h);
    }

    PotentialVals derivative2_r(double p, double r, double h = step) const {
        auto fr = [&](double x) { return this->operator()(p, x); };
        return stencil5second(fr, r, h);
    }

    virtual double Hamiltonian(double p, double q) const {
        auto vals = this->operator()(p, q);
        auto mass_eff =
            achilles::Constant::mN + vals.rscalar + std::complex<double>(0, 1) * vals.iscalar;
        return sqrt(p * p + pow(mass_eff, 2)).real() + vals.rvector;
    }

    double BindingMomentum(double r) const {
        auto func = [&](double p) { return Hamiltonian(p, r) - Constant::mN; };
        Brent brent(func);
        return brent.CalcRoot(0, 1000);
    }

    double EnergySpectrum(double r, double p) const { return Hamiltonian(p, r) - Constant::mN; }

    bool IsCaptured(double r, double mom) const { return Hamiltonian(mom, r) <= Constant::mN; }

    virtual bool IsRelativistic() const { return false; }

    // NOTE: Calculates in the non-relativistic limit
    // TODO: Modify to include scalar potential
    double Mstar(double p, double m, double r) const {
        return p / (p / m + derivative_p(p, r).rvector);
    }

    // NOTE: Calculates in the non-relativistic limit
    // TODO: Modify to include scalar potential
    double InMediumCorrectionNonRel(FourVector p1, FourVector p2, double m, double r1,
                                    double r2 = -1, double r3 = -1) const {
        if(r2 < 0) r2 = r1;
        if(r3 < 0) r3 = r1;
        double m1Star = Mstar(p1.P(), m, r1);
        double m2Star = Mstar(p2.P(), m, r2);
        double p12 = sqrt((p1.P2() + p2.P2()) / 2);
        double m12Star = Mstar(p12, m, r3);

        return (p1 - p2).P() / m / (p1 / m1Star - p2 / m2Star).P() * m12Star / m;
    }

  protected:
    PotentialVals stencil5(std::function<achilles::PotentialVals(double)> f, double x,
                           double h) const;
    PotentialVals stencil5second(std::function<achilles::PotentialVals(double)> f, double x,
                                 double h) const;
    std::array<PotentialVals, 3> stencil5all(std::function<achilles::PotentialVals(double)> f,
                                             double x, double h) const;

    static constexpr double step = 0.01;
};

template <typename Derived>
using RegistrablePotential =
    Registrable<Potential, Derived, std::shared_ptr<Nucleus> &, const YAML::Node &>;
using PotentialFactory = Factory<Potential, std::shared_ptr<Nucleus> &, const YAML::Node &>;

/*class SquareWellPotential : public Potential, RegistrablePotential<SquareWellPotential> {
    public:
        SquareWellPotential(std::shared_ptr<Nucleus> nucleus) : m_nucleus{std::move(nucleus)} {}

        PotentialVals operator()(const double&, const double &r) const override;
        static std::unique_ptr<Potential> Construct(std::shared_ptr<Nucleus>&, const YAML::Node&);
        static std::string Name() { return "SquareWellPotential"; }
        std::string GetReference() const override { return ""; }

    private:
        std::shared_ptr<Nucleus> m_nucleus;
};
*/

class WiringaPotential : public Potential, RegistrablePotential<WiringaPotential> {
  public:
    WiringaPotential(std::shared_ptr<Nucleus> nucleus, const double &rho0 = 0.16)
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

    static std::string Name() { return "Wiringa"; }
    static std::unique_ptr<Potential> Construct(std::shared_ptr<Nucleus> &, const YAML::Node &);

    double Hamiltonian(double p, double q) const override {
        auto vals = this->operator()(p, q);
        return Constant::mN + p * p / (2 * Constant::mN) + vals.rvector;
    }

    std::string GetReference() const override { return m_ref.GetReference(); }
    double Rho0() const { return m_rho0; }

    PotentialVals operator()(double plab, double radius) const override;

  private:
    std::shared_ptr<Nucleus> m_nucleus;
    double m_rho0;
    Reference m_ref;
};

class CooperPotential : public Potential, RegistrablePotential<CooperPotential> {
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
                data[8 * j + i] = pt1[j * 8 + i];
                data[8 * (j + 8) + i] = pt2[j * 8 + i];
                if(j < 6) data[8 * (j + 16) + i] = pt3[i * 6 + j];
            }
        }
    }

    static std::string Name() { return "Cooper"; }
    static std::unique_ptr<Potential> Construct(std::shared_ptr<Nucleus> &, const YAML::Node &);

    std::string GetReference() const override { return m_ref.GetReference(); }
    bool IsRelativistic() const override { return true; }

    PotentialVals operator()(double plab, double radius) const override {
        return evaluate(plab, radius);
    }

    PotentialVals evaluate(const double &plab, const double &radius) const;

  protected:
    std::shared_ptr<Nucleus> m_nucleus;

  private:
    Reference m_ref;
    std::array<double, 22 * 8> data{};

    double Data(size_t i, size_t j) const { return data[8 * (i - 1) + j - 1]; }
    double CalcTerm(double prefact, double real, double acb, double imag, double radius) const {
        const auto s1 = sech(real * acb / imag);
        const auto s2 = sech(radius / imag);
        const auto prod = s1 * s2;
        const auto a = s1 - prod;
        const auto b = s2 - prod;
        return prefact * b / (a + b);
    }
    double CalcTermSurf(double prefact, double real, double acb, double imag, double radius) const {
        const auto s1 = sech(real * acb / imag);
        const auto s2 = sech(radius / imag);
        const auto prod = s1 * s2;
        const auto a = s1 - prod;
        const auto b = s2 - prod;
        return prefact * a * b / (a + b) / (a + b);
    }

    static constexpr double wp = 1.0072545 * achilles::Constant::AMU;
    static constexpr std::array<double, 8 * 8> pt1{
        2.313828E+01,  -1.233380E+02, 2.432987E+02,  -2.092616E+02, 6.606549E+01,  5.645371E+00,
        -1.014396E+01, 5.733192E+00,  -6.841731E+00, 3.747771E+01,  -6.858889E+01, 5.484194E+01,
        -1.617735E+01, 9.973483E-01,  -1.044247E+00, 4.707822E-01,  1.688843E+01,  -8.404292E+01,
        1.550125E+02,  -1.251468E+02, 3.722361E+01,  2.912163E+00,  -3.560279E+00, 1.645079E+00,
        2.143528E+02,  -9.614898E+02, 1.576657E+03,  -1.147987E+03, 3.136021E+02,  2.226300E+01,
        -3.363301E+01, 1.704743E+01,  7.852790E+00,  -4.551995E+01, 9.525029E+01,  -8.644521E+01,
        2.891801E+01,  5.105014E+00,  -7.512490E+00, 3.681682E+00,  -8.332038E+00, 4.210331E+01,
        -7.878047E+01, 6.498050E+01,  -1.997590E+01, 4.603199E+00,  -7.773754E+00, 4.281708E+00,
        1.771733E+01,  -9.514249E+01, 1.826431E+02,  -1.507498E+02, 4.563918E+01,  5.566540E+00,
        -1.001441E+01, 5.551678E+00,  -6.557412E+00, 3.596972E+01,  -6.571700E+01, 5.233473E+01,
        -1.534723E+01, 1.007355E+00,  -9.669271E-01, 4.066760E-01};
    static constexpr std::array<double, 8 * 8> pt2{
        2.183389E+01,  -1.085667E+02, 2.012175E+02,  -1.637375E+02, 4.926579E+01,  2.874114E+00,
        -3.559180E+00, 1.647730E+00,  -5.944315E+01, 3.517880E+02,  -6.896441E+02, 6.080210E+02,
        -2.026092E+02, -3.492561E+01, 5.197784E+01,  -2.516421E+01, 6.414280E-01,  -9.685444E+00,
        2.992500E+01,  -3.446931E+01, 1.361126E+01,  5.209321E+00,  -7.563083E+00, 3.671130E+00,
        -7.989428E+00, 4.318908E+01,  -8.862678E+01, 7.860778E+01,  -2.582041E+01, 7.390220E+00,
        -1.207679E+01, 6.360783E+00,  1.871682E+01,  -1.144552E+02, 2.566304E+02,  -2.509193E+02,
        9.003415E+01,  0.000000E+00,  0.000000E+00,  0.000000E+00,  1.000000E+00,  0.000000E+00,
        0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
        1.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
        0.000000E+00,  0.000000E+00,  3.713738E+01,  -1.952369E+02, 3.842865E+02,  -3.379853E+02,
        1.116151E+02,  0.000000E+00,  0.000000E+00,  0.000000E+00};
    static constexpr std::array<double, 8 * 6> pt3{
        1.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00,
        0.000000E+00, 0.000000E+00, 1.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00,
        0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00,
        0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00,
        0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00,
        0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00,
        0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00,
        0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00};
};

class SchroedingerPotential : public CooperPotential, RegistrablePotential<SchroedingerPotential> {
  public:
    SchroedingerPotential(std::shared_ptr<Nucleus> nucleus, size_t mode)
        : CooperPotential(nucleus), m_mode{mode} {}

    static std::string Name() { return "Schroedinger"; }
    static std::unique_ptr<Potential> Construct(std::shared_ptr<Nucleus> &, const YAML::Node &);

    PotentialVals operator()(double plab, double radius) const override;

    double Hamiltonian(double p, double q) const override {
        auto vals = this->operator()(p, q);
        return Constant::mN + p * p / (2 * Constant::mN) + vals.rvector;
    }

  private:
    static constexpr double wp = 1.0072545 * achilles::Constant::AMU;
    static constexpr double hc2 = achilles::Constant::HBARC * achilles::Constant::HBARC;
    size_t m_mode;
};

} // namespace achilles

#endif
