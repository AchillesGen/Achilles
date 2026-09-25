// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef THREEVECTOR_HH
#define THREEVECTOR_HH

#include <array>
#include <cmath>
#include <iosfwd>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "Achilles/PhysicalUnits.hh"
#include "Achilles/UnitsFormat.hh"
#include "spdlog/fmt/ostr.h"

namespace achilles {

/// @brief ThreeVectorT is a container for dealing with three vectors
///
/// Templated on the mass dimension D of its components, so that one body serves
/// a three-momentum (D = 1, MeV), a position (D = -1, MeV^-1 -- ask for
/// .in(units::fm)) and a velocity (D = 0). Mixing them is a compile error.
template <int D> class ThreeVectorT {
  public:
    using RotMat = std::array<double, 9>;
    /// Component type: MeV^D
    using Q = units::Quantity<D>;
    /// Type of a product of two components: MeV^2D
    using Q2 = units::Quantity<2 * D>;

    /// @name Constructors and Destructors
    ///@{

    /// Create an empty ThreeVectorT object
    constexpr ThreeVectorT() noexcept : vec({Q{}, Q{}, Q{}}) {}
    /// Create a ThreeVectorT object with values given by p
    ///@param p: A std::array<Q, 3> containing the values for the vector
    constexpr explicit ThreeVectorT(std::array<Q, 3> p) noexcept : vec(p) {}
    /// Create a ThreeVectorT object with values pX, pY, and pZ
    ///@param pX: The x value of the vector
    ///@param pY: The y value of the vector
    ///@param pZ: The z value of the vector
    constexpr ThreeVectorT(Q pX, Q pY, Q pZ) noexcept : vec({pX, pY, pZ}) {}
    ThreeVectorT(const ThreeVectorT &other) noexcept = default;
    ThreeVectorT(ThreeVectorT &&other) noexcept = default;

    ThreeVectorT &operator=(const ThreeVectorT &) noexcept = default;
    ThreeVectorT &operator=(ThreeVectorT &&) noexcept = default;

    ~ThreeVectorT() = default;
    ///@}

    /// @name Interop boundary
    /// @{
    /// Raw canonical components (MeV^D). For a position this is MeV^-1, NOT fm;
    /// use ToArray(units::fm) at any fm-facing boundary.
    constexpr const std::array<Q, 3> &Native() const noexcept { return vec; }
    constexpr std::array<Q, 3> &Native() noexcept { return vec; }

    /// Build from raw canonical components (MeV^D)
    static constexpr ThreeVectorT FromNative(const std::array<double, 3> &p) noexcept {
        return ThreeVectorT{Q{p[0]}, Q{p[1]}, Q{p[2]}};
    }

    /// Components expressed in the given unit
    std::array<double, 3> ToArray(units::Unit<D> u) const noexcept {
        return {vec[0].in(u), vec[1].in(u), vec[2].in(u)};
    }
    /// @}

    /// @name Setters
    /// @{

    void SetXYZ(const std::array<Q, 3> &position) noexcept { vec = position; }
    void SetXYZ(Q x, Q y, Q z) noexcept { vec = std::array<Q, 3>{x, y, z}; }
    void SetPxPyPz(const std::array<Q, 3> &momentum) noexcept { vec = momentum; }
    void SetPxPyPz(Q pX, Q pY, Q pZ) noexcept { vec = std::array<Q, 3>{pX, pY, pZ}; }

    void SetX(Q x) noexcept { vec[0] = x; }
    void SetY(Q y) noexcept { vec[1] = y; }
    void SetZ(Q z) noexcept { vec[2] = z; }
    void SetPx(Q pX) noexcept { vec[0] = pX; }
    void SetPy(Q pY) noexcept { vec[1] = pY; }
    void SetPz(Q pZ) noexcept { vec[2] = pZ; }
    ///@}

    /// @name Getters
    /// @{

    constexpr size_t Size() const noexcept { return 3; }

    constexpr Q X() const noexcept { return vec[0]; }
    constexpr Q Y() const noexcept { return vec[1]; }
    constexpr Q Z() const noexcept { return vec[2]; }
    constexpr Q Px() const noexcept { return vec[0]; }
    constexpr Q Py() const noexcept { return vec[1]; }
    constexpr Q Pz() const noexcept { return vec[2]; }

    /// Transverse component squared
    constexpr Q2 Pt2() const noexcept { return vec[0] * vec[0] + vec[1] * vec[1]; }
    /// Transverse component
    Q Pt() const noexcept { return Q{std::sqrt(Pt2().native())}; }
    /// Magnitude squared
    constexpr Q2 P2() const noexcept { return (*this) * (*this); }
    /// Magnitude
    Q P() const noexcept { return Q{std::sqrt(P2().native())}; }
    constexpr Q2 Magnitude2() const noexcept { return P2(); }
    Q Magnitude() const noexcept { return P(); }

    /// Angle between the z-component and the transverse component (radians)
    double Theta() const noexcept { return std::atan2(Pt().native(), Pz().native()); }

    /// Angle between the x and y components (radians)
    double Phi() const noexcept {
        const double phi = std::atan2(Py().native(), Px().native());
        if(phi < 0) return phi + 2 * M_PI;
        return phi;
    }

    /// Rotate by three Euler angles
    ThreeVectorT Rotate(const std::array<double, 3> &angles) const noexcept {
        const double c1 = cos(angles[0]), s1 = sin(angles[0]);
        const double c2 = cos(angles[1]), s2 = sin(angles[1]);
        const double c3 = cos(angles[2]), s3 = sin(angles[2]);

        return {(c1 * c3 - c2 * s1 * s3) * vec[0] + (-c1 * s3 - c2 * c3 * s1) * vec[1] +
                    s1 * s2 * vec[2],
                (c3 * s1 + c1 * c2 * s3) * vec[0] + (c1 * c2 * c3 - s1 * s3) * vec[1] -
                    c1 * s2 * vec[2],
                s2 * s3 * vec[0] + c3 * s2 * vec[1] + c2 * vec[2]};
    }

    /// Rotate by a rotation matrix
    ThreeVectorT Rotate(const RotMat &mat) const noexcept {
        return {mat[0] * vec[0] + mat[1] * vec[1] + mat[2] * vec[2],
                mat[3] * vec[0] + mat[4] * vec[1] + mat[5] * vec[2],
                mat[6] * vec[0] + mat[7] * vec[1] + mat[8] * vec[2]};
    }

    /// Rotation matrix aligning this vector with the given (unit) axis
    RotMat Align(const ThreeVectorT<0> &axis) const noexcept {
        ThreeVectorT<0> a = Unit();

        auto v = a.Cross(axis);
        double c = a.Dot(axis).native();

        return {
            1 - v[1].native() * v[1].native() / (1 + c) - v[2].native() * v[2].native() / (1 + c),
            -v[2].native() + v[0].native() * v[1].native() / (1 + c),
            v[1].native() + v[0].native() * v[2].native() / (1 + c),
            v[2].native() + v[0].native() * v[1].native() / (1 + c),
            1 - v[0].native() * v[0].native() / (1 + c) - v[2].native() * v[2].native() / (1 + c),
            -v[0].native() + v[1].native() * v[2].native() / (1 + c),
            -v[1].native() + v[0].native() * v[2].native() / (1 + c),
            v[0].native() + v[1].native() * v[2].native() / (1 + c),
            1 - v[0].native() * v[0].native() / (1 + c) - v[1].native() * v[1].native() / (1 + c)};
    }

    /// Rotation matrix aligning this vector with the z-axis
    RotMat AlignZ() const noexcept {
        ThreeVectorT<0> z{units::Dimensionless{0}, units::Dimensionless{0},
                          units::Dimensionless{1}};
        return Align(z);
    }
    ///@}

    /// @name Functions
    /// @{

    /// Dot product; the dimensions add
    template <int B>
    constexpr units::Quantity<D + B> Dot(const ThreeVectorT<B> &other) const noexcept {
        return (*this) * other;
    }

    /// Cross product; the dimensions add
    template <int B>
    constexpr ThreeVectorT<D + B> Cross(const ThreeVectorT<B> &other) const noexcept {
        return {vec[1] * other[2] - vec[2] * other[1], vec[2] * other[0] - vec[0] * other[2],
                vec[0] * other[1] - vec[1] * other[0]};
    }

    /// Unit vector in the direction of this vector. Dimensionless by construction.
    ThreeVectorT<0> Unit() const noexcept {
        const Q norm = Magnitude();
        return {vec[0] / norm, vec[1] / norm, vec[2] / norm};
    }

    std::string ToString() const noexcept {
        return "ThreeVector(" + std::to_string(vec[0].native()) + ", " +
               std::to_string(vec[1].native()) + ", " + std::to_string(vec[2].native()) + ")";
    }
    ///@}

    /// @name Operator Overloads
    /// @{

    ThreeVectorT &operator+=(const ThreeVectorT &other) noexcept {
        vec[0] += other.vec[0];
        vec[1] += other.vec[1];
        vec[2] += other.vec[2];
        return *this;
    }

    ThreeVectorT &operator-=(const ThreeVectorT &other) noexcept {
        vec[0] -= other.vec[0];
        vec[1] -= other.vec[1];
        vec[2] -= other.vec[2];
        return *this;
    }

    ThreeVectorT &operator*=(double scale) noexcept {
        vec[0] *= scale;
        vec[1] *= scale;
        vec[2] *= scale;
        return *this;
    }

    ThreeVectorT &operator/=(double scale) {
        vec[0] /= scale;
        vec[1] /= scale;
        vec[2] /= scale;
        return *this;
    }

    /// Dot product; the dimensions add
    template <int B>
    constexpr units::Quantity<D + B> operator*(const ThreeVectorT<B> &other) const noexcept {
        return vec[0] * other[0] + vec[1] * other[1] + vec[2] * other[2];
    }

    ThreeVectorT operator-() const noexcept { return {-vec[0], -vec[1], -vec[2]}; }
    ThreeVectorT operator+() const noexcept { return *this; }

    ThreeVectorT operator*(double scale) const noexcept { return ThreeVectorT(*this) *= scale; }
    ThreeVectorT operator/(double scale) const { return ThreeVectorT(*this) /= scale; }

    /// Scale by a quantity; the dimensions add
    template <int B>
    constexpr ThreeVectorT<D + B> operator*(units::Quantity<B> scale) const noexcept {
        return {vec[0] * scale, vec[1] * scale, vec[2] * scale};
    }

    /// Divide by a quantity; the dimensions subtract
    template <int B>
    constexpr ThreeVectorT<D - B> operator/(units::Quantity<B> scale) const noexcept {
        return {vec[0] / scale, vec[1] / scale, vec[2] / scale};
    }

    ThreeVectorT operator-(const ThreeVectorT &other) const noexcept {
        return ThreeVectorT(*this) -= other;
    }
    ThreeVectorT operator+(const ThreeVectorT &other) const noexcept {
        return ThreeVectorT(*this) += other;
    }

    bool operator==(const ThreeVectorT &other) const noexcept { return vec == other.vec; }
    bool operator!=(const ThreeVectorT &other) const noexcept { return !(*this == other); }

    constexpr Q &operator[](const std::size_t &idx) { return vec[idx]; }
    constexpr const Q &operator[](const std::size_t &idx) const { return vec[idx]; }

    Q &at(const std::size_t &idx) {
        if(idx > 2) throw std::range_error("Max value is 2.");
        return vec[idx];
    }
    const Q &at(const std::size_t &idx) const {
        if(idx > 2) throw std::range_error("Max value is 2.");
        return vec[idx];
    }
    /// @}

    ///@name Stream Operators
    /// @{
    template <typename OStream> friend OStream &operator<<(OStream &os, const ThreeVectorT &vec3) {
        os << "ThreeVector(" << vec3.Px().native() << ", " << vec3.Py().native() << ", "
           << vec3.Pz().native() << ")";
        return os;
    }
    /// @}

  private:
    std::array<Q, 3> vec;
};

/// The storage is bit-compatible with the legacy double[3] layout, so ROOT,
/// HepMC3 and Fortran buffers are unaffected.
static_assert(sizeof(ThreeVectorT<1>) == 3 * sizeof(double));
static_assert(std::is_trivially_copyable_v<ThreeVectorT<1>>);
static_assert(std::is_standard_layout_v<ThreeVectorT<1>>);

template <int D> ThreeVectorT<D> operator*(double s, const ThreeVectorT<D> &v) noexcept {
    return v * s;
}

template <int B, int D>
constexpr ThreeVectorT<B + D> operator*(units::Quantity<B> s, const ThreeVectorT<D> &v) noexcept {
    return v * s;
}

template <int D> std::istream &operator>>(std::istream &is, ThreeVectorT<D> &vec) {
    std::string head(12, ' '), sep1(1, ' '), sep2(1, ' '), tail(1, ' ');
    double px{}, py{}, pz{};
    is.read(&head[0], 12);
    is >> px;
    is.read(&sep1[0], 1);
    is >> py;
    is.read(&sep2[0], 1);
    is >> pz;
    is.read(&tail[0], 1);
    if(head == "ThreeVector(" && sep1 == "," && sep2 == "," && tail == ")")
        vec = ThreeVectorT<D>::FromNative({px, py, pz});
    return is;
}

/// Three-momentum, canonical MeV
using ThreeVector = ThreeVectorT<1>;
/// Position, canonical MeV^-1 -- .in(units::fm) at any fm-facing boundary
using ThreePosition = ThreeVectorT<-1>;
/// Velocity / boost vector, dimensionless
using ThreeBoost = ThreeVectorT<0>;

} // namespace achilles

template <int D> struct fmt::formatter<achilles::ThreeVectorT<D>> {
    char presentation = 'e';
    constexpr auto parse(format_parse_context &ctx) -> format_parse_context::iterator {
        auto it = ctx.begin(), end = ctx.end();
        if(it != end && (*it == 'f' || *it == 'e')) presentation = *it++;
        if(it != end && *it != '}') throw format_error("Invalid format");
        return it;
    }

    auto format(const achilles::ThreeVectorT<D> &p, format_context &ctx) const
        -> format_context::iterator {
        return format_to(ctx.out(),
                         presentation == 'f' ? "ThreeVector({:.8f}, {:.8f}, {:.8f})"
                                             : "ThreeVector({:.8e}, {:.8e}, {:.8e})",
                         p.Px().native(), p.Py().native(), p.Pz().native());
    }
};

#endif // end of include guard: THREEVECTOR_HH
