// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef FOURVECTOR_HH
#define FOURVECTOR_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <iosfwd>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "Achilles/PhysicalUnits.hh"
#include "Achilles/ThreeVector.hh"
#include "fmt/format.h"
#include "spdlog/fmt/ostr.h"

namespace achilles {

/// @brief FourVectorT is a container for dealing with four vectors
///
/// Templated on the mass dimension D of its components, so that one body serves
/// both a four-momentum (D = 1, MeV) and a four-position (D = -1, MeV^-1 --
/// ask for .in(units::fm)). Adding one to the other is a compile error.
template <int D> class FourVectorT {
  public:
    using RotMat = std::array<double, 9>;
    /// Component type: MeV^D
    using Q = units::Quantity<D>;
    /// Type of a product of two components: MeV^2D
    using Q2 = units::Quantity<2 * D>;

    /// @name Constructors and Destructors
    ///@{

    /// Create an empty FourVectorT object
    constexpr FourVectorT() noexcept : vec({Q{}, Q{}, Q{}, Q{}}) {}
    /// Create a FourVectorT object with values given by p
    ///@param p: A std::array<Q, 4> containing the values for the vector
    constexpr explicit FourVectorT(std::array<Q, 4> p) noexcept : vec(p) {}
    /// Create a FourVectorT object with values E, pX, pY, and pZ
    constexpr FourVectorT(Q E, Q pX, Q pY, Q pZ) noexcept : vec({E, pX, pY, pZ}) {}
    /// Create a FourVectorT object from a ThreeVectorT and an energy
    constexpr FourVectorT(const ThreeVectorT<D> &other, Q E) noexcept
        : vec({E, other[0], other[1], other[2]}) {}
    FourVectorT(const FourVectorT &other) noexcept = default;
    FourVectorT(FourVectorT &&) noexcept = default;

    ~FourVectorT() = default;
    ///@}

    /// @name Interop boundary
    /// @{
    /// Raw canonical components (MeV^D), ordered {E, px, py, pz}. For a position
    /// this is MeV^-1, NOT fm; use ToArray(units::fm) at any fm-facing boundary.
    constexpr const std::array<Q, 4> &Native() const noexcept { return vec; }
    constexpr std::array<Q, 4> &Native() noexcept { return vec; }

    /// Build from raw canonical components (MeV^D)
    static constexpr FourVectorT FromNative(const std::array<double, 4> &p) noexcept {
        return FourVectorT{Q{p[0]}, Q{p[1]}, Q{p[2]}, Q{p[3]}};
    }

    /// Components expressed in the given unit
    std::array<double, 4> ToArray(units::Unit<D> u) const noexcept {
        return {vec[0].in(u), vec[1].in(u), vec[2].in(u), vec[3].in(u)};
    }
    /// @}

    /// @name Setters
    /// @{

    void SetPxPyPzE(const std::array<Q, 4> &p) noexcept { vec = p; }
    void SetPxPyPzE(Q pX, Q pY, Q pZ, Q E) noexcept { vec = std::array<Q, 4>{E, pX, pY, pZ}; }

    /// Set the four vector given a three vector and an invariant mass
    void SetVectM(const ThreeVectorT<D> &p, Q mass) noexcept {
        Px() = p.Px();
        Py() = p.Py();
        Pz() = p.Pz();
        E() = Q{std::sqrt((mass * mass + p * p).native())};
    }

    void SetPx(Q pX) noexcept { vec[1] = pX; }
    void SetPy(Q pY) noexcept { vec[2] = pY; }
    void SetPz(Q pZ) noexcept { vec[3] = pZ; }
    void SetE(Q E) noexcept { vec[0] = E; }
    ///@}

    /// @name Getters
    /// @{

    constexpr size_t Size() const noexcept { return 4; }

    constexpr Q X() const noexcept { return vec[1]; }
    constexpr Q Y() const noexcept { return vec[2]; }
    constexpr Q Z() const noexcept { return vec[3]; }
    constexpr Q T() const noexcept { return vec[0]; }

    constexpr Q Px() const noexcept { return vec[1]; }
    constexpr Q &Px() noexcept { return vec[1]; }
    constexpr Q Py() const noexcept { return vec[2]; }
    constexpr Q &Py() noexcept { return vec[2]; }
    constexpr Q Pz() const noexcept { return vec[3]; }
    constexpr Q &Pz() noexcept { return vec[3]; }
    constexpr Q E() const noexcept { return vec[0]; }
    constexpr Q &E() noexcept { return vec[0]; }

    /// Transverse momentum squared
    constexpr Q2 Pt2() const noexcept { return Px() * Px() + Py() * Py(); }
    /// Transverse momentum
    Q Pt() const noexcept { return Q{std::sqrt(Pt2().native())}; }
    /// Three momentum squared
    constexpr Q2 P2() const noexcept { return Pt2() + Pz() * Pz(); }
    /// Three momentum
    Q P() const noexcept { return Q{std::sqrt(P2().native())}; }

    /// Invariant mass squared
    constexpr Q2 M2() const noexcept { return (*this) * (*this); }

    /// Invariant mass
    Q M() const noexcept {
        const double mass2 = M2().native();
        if(std::abs(mass2) < tolerance) return Q{0.0};
        return Q{std::sqrt(mass2)};
    }

    constexpr Q2 Magnitude2() const noexcept { return M2(); }
    Q Magnitude() const noexcept { return M(); }

    /// Angle between the z-axis and the transverse plane (radians)
    double Theta() const noexcept { return std::atan2(Pt().native(), Pz().native()); }
    double CosTheta() const noexcept { return std::cos(Theta()); }

    /// Angle in the transverse plane (radians)
    double Phi() const noexcept {
        const double phi = std::atan2(Py().native(), Px().native());
        if(phi < 0) return phi + 2 * M_PI;
        return phi;
    }

    /// Rapidity (dimensionless)
    double Rapidity() const noexcept {
        return std::log(((E() + Pz()) / (E() - Pz())).native()) / 2;
    }

    /// Distance in the Eta-Phi plane between two four vectors
    double DeltaR(const FourVectorT &other) const noexcept {
        const double DEta = Rapidity() - other.Rapidity();
        const double DPhi = Phi() - other.Phi();
        return std::sqrt(DEta * DEta + DPhi * DPhi);
    }

    /// The spatial component
    ThreeVectorT<D> Vec3() const noexcept { return {Px(), Py(), Pz()}; }
    ///@}

    /// @name Functions
    /// @{

    /// 1 - cos(theta) evaluated without cancellation for small angles
    double SmallOMCT(const FourVectorT &v) const noexcept {
        const double mag = std::sqrt((P2() * v.P2()).native());
        const double pq = (vec[1] * v[1] + vec[2] * v[2] + vec[3] * v[3]).native();
        const double ct = std::min(std::max(pq / mag, -1.), 1.);
        if(ct < 0.) return 1. - ct;
        const double st = Vec3().Cross(v.Vec3()).P().native() / mag;
        const double st2 = st / (2. * std::sqrt((ct + 1) / 2.));
        return 2. * st2 * st2;
    }

    Q2 SmallMLDP(const FourVectorT &v) const noexcept { return vec[0] * v[0] * SmallOMCT(v); }

    /// Boost by a dimensionless velocity. Passing a momentum here will not compile.
    FourVectorT Boost(const ThreeVectorT<0> &beta) const noexcept {
        const double beta2 = (beta * beta).native();
        const double gamma = 1.0 / std::sqrt(1.0 - beta2);
        const Q betap = beta[0] * Px() + beta[1] * Py() + beta[2] * Pz();
        const double gamma2 = beta2 > 0 ? (gamma - 1.0) / beta2 : 0.0;

        return {gamma * (E() + betap), Px() + gamma2 * betap * beta[0] + gamma * beta[0] * E(),
                Py() + gamma2 * betap * beta[1] + gamma * beta[1] * E(),
                Pz() + gamma2 * betap * beta[2] + gamma * beta[2] * E()};
    }

    FourVectorT Boost(double beta_x, double beta_y, double beta_z) const noexcept {
        return Boost(ThreeVectorT<0>{units::Dimensionless{beta_x}, units::Dimensionless{beta_y},
                                     units::Dimensionless{beta_z}});
    }

    /// Rotate the spatial part by a rotation matrix
    FourVectorT Rotate(const RotMat &mat) const noexcept {
        return {E(), mat[0] * Px() + mat[1] * Py() + mat[2] * Pz(),
                mat[3] * Px() + mat[4] * Py() + mat[5] * Pz(),
                mat[6] * Px() + mat[7] * Py() + mat[8] * Pz()};
    }

    FourVectorT RotateBack(const RotMat &mat) const noexcept {
        return {E(), mat[0] * Px() + mat[3] * Py() + mat[6] * Pz(),
                mat[1] * Px() + mat[4] * Py() + mat[7] * Pz(),
                mat[2] * Px() + mat[5] * Py() + mat[8] * Pz()};
    }

    /// Rotation matrix aligning the spatial part with the given (unit) axis
    RotMat Align(const ThreeVectorT<0> &axis) const noexcept { return Vec3().Align(axis); }

    /// Rotation matrix aligning the spatial part with the z-axis
    RotMat AlignZ() const noexcept { return Vec3().AlignZ(); }

    /// Cross product of the spatial parts; the dimensions add
    template <int B> FourVectorT<D + B> Cross(const FourVectorT<B> &other) const noexcept {
        return {Vec3().Cross(other.Vec3()), units::Quantity<D + B>{0}};
    }

    /// The boost required to reach the rest frame of this vector. Dimensionless.
    ThreeVectorT<0> BoostVector() const noexcept { return {Px() / E(), Py() / E(), Pz() / E()}; }

    /// Minkowski dot product; the dimensions add
    template <int B>
    constexpr units::Quantity<D + B> Dot(const FourVectorT<B> &other) const noexcept {
        return (*this) * other;
    }

    /// Cosine of the angle between the spatial parts
    double CosAngle(const FourVectorT &other) const noexcept {
        auto p1 = Vec3(), p2 = other.Vec3();
        return std::max(std::min(((p1 * p2) / (p1.P() * p2.P())).native(), 1.0), -1.0);
    }

    double Angle(const FourVectorT &other) const noexcept { return std::acos(CosAngle(other)); }

    std::string ToString() const noexcept {
        return "FourVector(" + std::to_string(vec[0].native()) + ", " +
               std::to_string(vec[1].native()) + ", " + std::to_string(vec[2].native()) + ", " +
               std::to_string(vec[3].native()) + ")";
    }
    ///@}

    /// @name Operator Overloads
    /// @{

    FourVectorT &operator=(const FourVectorT &) noexcept = default;
    FourVectorT &operator=(FourVectorT &&) noexcept = default;

    FourVectorT &operator+=(const FourVectorT &other) noexcept {
        vec[0] += other.vec[0];
        vec[1] += other.vec[1];
        vec[2] += other.vec[2];
        vec[3] += other.vec[3];
        return *this;
    }

    FourVectorT &operator-=(const FourVectorT &other) noexcept {
        vec[0] -= other.vec[0];
        vec[1] -= other.vec[1];
        vec[2] -= other.vec[2];
        vec[3] -= other.vec[3];
        return *this;
    }

    FourVectorT &operator*=(double scale) noexcept {
        vec[0] *= scale;
        vec[1] *= scale;
        vec[2] *= scale;
        vec[3] *= scale;
        return *this;
    }

    FourVectorT &operator/=(double scale) {
        vec[0] /= scale;
        vec[1] /= scale;
        vec[2] /= scale;
        vec[3] /= scale;
        return *this;
    }

    /// Minkowski dot product; the dimensions add
    template <int B>
    constexpr units::Quantity<D + B> operator*(const FourVectorT<B> &other) const noexcept {
        return vec[0] * other[0] - (vec[1] * other[1] + vec[2] * other[2] + vec[3] * other[3]);
    }

    FourVectorT operator-() const noexcept { return {-vec[0], -vec[1], -vec[2], -vec[3]}; }
    FourVectorT operator+() const noexcept { return *this; }

    FourVectorT operator*(double scale) const noexcept { return FourVectorT(*this) *= scale; }
    FourVectorT operator/(double scale) const { return FourVectorT(*this) /= scale; }

    /// Scale by a quantity; the dimensions add
    template <int B>
    constexpr FourVectorT<D + B> operator*(units::Quantity<B> scale) const noexcept {
        return {vec[0] * scale, vec[1] * scale, vec[2] * scale, vec[3] * scale};
    }

    /// Divide by a quantity; the dimensions subtract
    template <int B>
    constexpr FourVectorT<D - B> operator/(units::Quantity<B> scale) const noexcept {
        return {vec[0] / scale, vec[1] / scale, vec[2] / scale, vec[3] / scale};
    }

    FourVectorT operator-(const FourVectorT &other) const noexcept {
        return FourVectorT(*this) -= other;
    }
    FourVectorT operator+(const FourVectorT &other) const noexcept {
        return FourVectorT(*this) += other;
    }

    bool operator==(const FourVectorT &other) const noexcept { return vec == other.vec; }
    bool operator!=(const FourVectorT &other) const noexcept { return !(*this == other); }

    /// Approximate equality, with the tolerance given in canonical units (MeV^D)
    bool Approx(const FourVectorT &other, double eps = 1e-8) const noexcept {
        for(size_t i = 0; i < vec.size(); ++i) {
            if(std::abs(vec[i].native() - other.vec[i].native()) > eps) return false;
        }
        return true;
    }

    constexpr Q &operator[](const std::size_t &idx) { return vec[idx]; }
    constexpr const Q &operator[](const std::size_t &idx) const { return vec[idx]; }

    Q &at(const std::size_t &idx) {
        if(idx > 3) throw std::range_error("Max value is 3.");
        return vec[idx];
    }
    const Q &at(const std::size_t &idx) const {
        if(idx > 3) throw std::range_error("Max value is 3.");
        return vec[idx];
    }
    /// @}

    ///@name Stream Operators
    /// @{
    template <typename OStream> friend OStream &operator<<(OStream &os, const FourVectorT &vec4) {
        os << "FourVector(" << vec4.E().native() << ", " << vec4.Px().native() << ", "
           << vec4.Py().native() << ", " << vec4.Pz().native() << ")";
        return os;
    }
    /// @}

  private:
    std::array<Q, 4> vec;
    static constexpr double tolerance = 1e-12;
};

/// The storage is bit-compatible with the legacy double[4] layout, so ROOT,
/// HepMC3 and Fortran buffers are unaffected.
static_assert(sizeof(FourVectorT<1>) == 4 * sizeof(double));
static_assert(std::is_trivially_copyable_v<FourVectorT<1>>);
static_assert(std::is_standard_layout_v<FourVectorT<1>>);

template <int D> FourVectorT<D> operator*(double s, const FourVectorT<D> &v) noexcept {
    return v * s;
}

template <int B, int D>
constexpr FourVectorT<B + D> operator*(units::Quantity<B> s, const FourVectorT<D> &v) noexcept {
    return v * s;
}

template <int D> std::istream &operator>>(std::istream &is, FourVectorT<D> &vec) {
    std::string head_name = "FourVector(";
    std::string head(head_name.size(), ' '), sep1(1, ' '), sep2(1, ' '), sep3(1, ' '), tail(1, ' ');
    double e{}, px{}, py{}, pz{};
    is.read(&head[0], 11);
    is >> e;
    is.read(&sep1[0], 1);
    is >> px;
    is.read(&sep2[0], 1);
    is >> py;
    is.read(&sep3[0], 1);
    is >> pz;
    is.read(&tail[0], 1);
    if(head == head_name && sep1 == "," && sep2 == "," && sep3 == "," && tail == ")")
        vec = FourVectorT<D>::FromNative({e, px, py, pz});
    return is;
}

/// Four-momentum, canonical MeV
using FourVector = FourVectorT<1>;
/// Four-position, canonical MeV^-1 -- .in(units::fm) at any fm-facing boundary
using FourPosition = FourVectorT<-1>;

} // namespace achilles

template <int D> struct fmt::formatter<achilles::FourVectorT<D>> {
    char presentation = 'e';
    constexpr auto parse(format_parse_context &ctx) -> format_parse_context::iterator {
        auto it = ctx.begin(), end = ctx.end();
        if(it != end && (*it == 'f' || *it == 'e')) presentation = *it++;
        if(it != end && *it != '}') throw format_error("Invalid format");
        return it;
    }

    auto format(const achilles::FourVectorT<D> &p, format_context &ctx) const
        -> format_context::iterator {
        return format_to(ctx.out(),
                         presentation == 'f' ? "FourVector({:.8f}, {:.8f}, {:.8f}, {:.8f})"
                                             : "FourVector({:.8e}, {:.8e}, {:.8e}, {:.8e})",
                         p.E().native(), p.Px().native(), p.Py().native(), p.Pz().native());
    }
};

#endif /* end of include guard: FOURVECTOR_HH */
