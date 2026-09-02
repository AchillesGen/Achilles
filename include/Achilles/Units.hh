// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Strong-typed natural units for Achilles (single mass-dimension axis).
//
// In natural units (hbar = c = 1) every dimensionful quantity is a power of
// energy. We store ONE canonical unit -- MeV -- and tag each quantity with its
// mass dimension D (the power of energy):
//
//     Energy / momentum / mass   D = +1   stored in MeV
//     Length / time / 1/energy   D = -1   stored in MeV^-1
//     Area / cross section       D = -2   stored in MeV^-2
//     Energy^2 (e.g. Mandelstam) D = +2   stored in MeV^2
//     Dimensionless              D =  0
//
// Consequences of the single axis (all physically correct in natural units):
//   * A length and a cross section are just Quantity<-1> and Quantity<-2>;
//     a geometric area (fm^2) and a matrix-element cross section (MeV^-2) are
//     the SAME type and interconvert with no explicit hbar*c in user code.
//   * hbar*c = 197.3269804 MeV*fm is the ONE bridge constant. It is not a
//     dimensionful multiplier here -- it is dimensionless (= 1), so it appears
//     only inside the fm/barn unit *scales*, defined once below.
//   * A Quantity<-1> cannot know whether it is a length or a lifetime, so it
//     stores and prints as MeV^-1; ask for .in(fm) or .in(s) to choose.

#ifndef UNITS_HH
#define UNITS_HH

#include <cmath>
#include <string_view>
#include <type_traits>

namespace achilles::units {

template <int D> class Unit; // defined below

// ---------------------------------------------------------------------------
// Quantity<D>: a double stored in canonical units (MeV^D), tagged by mass
// dimension D. Zero runtime overhead; D exists only at compile time.
// ---------------------------------------------------------------------------
template <int D> class Quantity {
    double val_{};

  public:
    static constexpr int dim = D;

    constexpr Quantity() = default;
    explicit constexpr Quantity(double v) noexcept : val_(v) {}

    // Canonical value (MeV^D). Prefer .in(unit) at ABI boundaries -- for a
    // length this is MeV^-1, NOT fm.
    constexpr double native() const noexcept { return val_; }

    // Same-dimension add / sub / negate --------------------------------------
    friend constexpr Quantity operator+(Quantity a, Quantity b) noexcept {
        return Quantity{a.val_ + b.val_};
    }
    friend constexpr Quantity operator-(Quantity a, Quantity b) noexcept {
        return Quantity{a.val_ - b.val_};
    }
    constexpr Quantity operator-() const noexcept { return Quantity{-val_}; }
    constexpr Quantity &operator+=(Quantity o) noexcept {
        val_ += o.val_;
        return *this;
    }
    constexpr Quantity &operator-=(Quantity o) noexcept {
        val_ -= o.val_;
        return *this;
    }

    // Scale by a plain number -------------------------------------------------
    friend constexpr Quantity operator*(Quantity q, double s) noexcept {
        return Quantity{q.val_ * s};
    }
    friend constexpr Quantity operator*(double s, Quantity q) noexcept {
        return Quantity{s * q.val_};
    }
    friend constexpr Quantity operator/(Quantity q, double s) noexcept {
        return Quantity{q.val_ / s};
    }
    constexpr Quantity &operator*=(double s) noexcept {
        val_ *= s;
        return *this;
    }
    constexpr Quantity &operator/=(double s) noexcept {
        val_ /= s;
        return *this;
    }

    // Same-dimension comparisons ---------------------------------------------
    friend constexpr bool operator<(Quantity a, Quantity b) noexcept { return a.val_ < b.val_; }
    friend constexpr bool operator>(Quantity a, Quantity b) noexcept { return a.val_ > b.val_; }
    friend constexpr bool operator<=(Quantity a, Quantity b) noexcept { return a.val_ <= b.val_; }
    friend constexpr bool operator>=(Quantity a, Quantity b) noexcept { return a.val_ >= b.val_; }
    friend constexpr bool operator==(Quantity a, Quantity b) noexcept { return a.val_ == b.val_; }
    friend constexpr bool operator!=(Quantity a, Quantity b) noexcept { return a.val_ != b.val_; }

    // Extract in a named unit of the SAME dimension --------------------------
    template <int D2> constexpr double in(Unit<D2> u) const noexcept;
};

// The wrapper must stay layout-identical to the double it replaces so that
// arrays of Quantity can alias the buffers handed to Fortran/ROOT/HepMC3.
static_assert(sizeof(Quantity<1>) == sizeof(double), "Quantity must be a plain double");
static_assert(std::is_trivially_copyable<Quantity<1>>::value, "Quantity must be memcpy-able");
static_assert(std::is_standard_layout<Quantity<1>>::value, "Quantity must be standard layout");

// Dimension-combining multiply / divide -------------------------------------
template <int A, int B> constexpr Quantity<A + B> operator*(Quantity<A> a, Quantity<B> b) noexcept {
    return Quantity<A + B>{a.native() * b.native()};
}
template <int A, int B> constexpr Quantity<A - B> operator/(Quantity<A> a, Quantity<B> b) noexcept {
    return Quantity<A - B>{a.native() / b.native()};
}
template <int D> constexpr Quantity<-D> operator/(double s, Quantity<D> q) noexcept {
    return Quantity<-D>{s / q.native()};
}

namespace detail {
// constexpr-friendly Newton iteration, used only during constant evaluation.
constexpr double sqrt_constexpr(double x) noexcept {
    if(!(x > 0)) return x; // 0, -0 and NaN pass through; negatives handled below
    double r = x;
    for(int i = 0; i < 60; ++i) r = 0.5 * (r + x / r);
    return r;
}
constexpr double sqrt_impl(double x) noexcept {
#if defined(__has_builtin)
#if __has_builtin(__builtin_is_constant_evaluated)
    if(!__builtin_is_constant_evaluated()) return std::sqrt(x);
#endif
#endif
    return sqrt_constexpr(x);
}
} // namespace detail

// sqrt halves the dimension (even exponents only) ---------------------------
template <int D> constexpr Quantity<D / 2> sqrt(Quantity<D> q) noexcept {
    static_assert(D % 2 == 0, "sqrt of a Quantity requires an even mass dimension");
    return Quantity<D / 2>{detail::sqrt_impl(q.native())};
}

template <int D> constexpr Quantity<D> abs(Quantity<D> q) noexcept {
    return q.native() < 0 ? -q : q;
}

// ---------------------------------------------------------------------------
// Unit<D>: a labelled scale factor. One 'unit' equals `scale` canonical MeV^D.
// ---------------------------------------------------------------------------
template <int D> class Unit {
    double scale_;
    std::string_view name_;

  public:
    explicit constexpr Unit(double s, std::string_view n = {}) noexcept : scale_(s), name_(n) {}
    constexpr double scale() const noexcept { return scale_; }
    constexpr std::string_view name() const noexcept { return name_; }
};

template <int D> constexpr Quantity<D> operator*(double m, Unit<D> u) noexcept {
    return Quantity<D>{m * u.scale()};
}
template <int D> constexpr Quantity<D> operator*(Unit<D> u, double m) noexcept {
    return Quantity<D>{m * u.scale()};
}

template <int D> template <int D2> constexpr double Quantity<D>::in(Unit<D2> u) const noexcept {
    static_assert(D == D2, "unit mismatch: this quantity cannot be expressed in that unit");
    return val_ / u.scale();
}

// ---------------------------------------------------------------------------
// Named dimensions
// ---------------------------------------------------------------------------
using Dimensionless = Quantity<0>;
using Energy = Quantity<1>; // == Momentum == Mass  (c = 1)
using Momentum = Quantity<1>;
using Mass = Quantity<1>;
using Energy2 = Quantity<2>;
using InvEnergy = Quantity<-1>; // == Length == Time    (hbar = c = 1)
using Length = Quantity<-1>;
using Time = Quantity<-1>;
using Area = Quantity<-2>; // == CrossSection
using CrossSection = Quantity<-2>;

// ---------------------------------------------------------------------------
// The single bridge constant. Everything on the length/area axes derives from
// it, so fm, barn, and HBARC2 can never fall out of sync.
//   hbar*c = 197.32698045930246 MeV*fm = 1  =>  1 fm = (1/hbar c) MeV^-1
//
// hbar*c is not written down as a literal: it is built from the SI-exact 2019
// values of c and h exactly as Constants.hh has always built it, so every
// number the generator produces is unchanged. Constant::C / H / HBAR / HBARC /
// HBARC2 are aliases of these.
// ---------------------------------------------------------------------------
inline constexpr double kSpeedOfLight = 2.99792458e8 * 1.0e15;      // fm / s (exact)
inline constexpr double kPlanck = 6.62607015e-34 / 1.602176634e-13; // MeV s  (exact)
inline constexpr double kHbar = kPlanck / (2 * M_PI);               // MeV s
inline constexpr double kHbarC = kHbar * kSpeedOfLight;             // MeV fm
inline constexpr double kHbarC2_MeV2mb = kHbarC * kHbarC * 10;      // mb MeV^2
inline constexpr double kFm = 1.0 / kHbarC;                         // MeV^-1 per fm
inline constexpr double kBarn = 100.0 * kFm * kFm; // MeV^-2 per barn (1 b = 100 fm^2)

// Energy units (D = +1), canonical MeV -- the MeV/GeV fix.
inline constexpr Unit<1> eV{1.0e-6, "eV"};
inline constexpr Unit<1> keV{1.0e-3, "keV"};
inline constexpr Unit<1> MeV{1.0, "MeV"};
inline constexpr Unit<1> GeV{1.0e3, "GeV"};
inline constexpr Unit<1> TeV{1.0e6, "TeV"};
inline constexpr Unit<2> MeV2{1.0, "MeV^2"};
inline constexpr Unit<2> GeV2{1.0e6, "GeV^2"};

// Length / time units (D = -1), canonical MeV^-1. hbar*c lives in the scale.
inline constexpr Unit<-1> iMeV{1.0, "MeV^-1"};
inline constexpr Unit<-1> iGeV{1.0e-3, "GeV^-1"}; // 1 GeV^-1 = 1e-3 MeV^-1
inline constexpr Unit<-1> fm{kFm, "fm"};
inline constexpr Unit<-1> nm{kFm * 1.0e6, "nm"};
inline constexpr Unit<-1> m{kFm * 1.0e15, "m"};

// Area / cross-section units (D = -2), canonical MeV^-2. Barn family derives
// from the fm scale, so MeV^-2 <-> nb <-> fm^2 are automatically consistent.
inline constexpr Unit<-2> iMeV2{1.0, "MeV^-2"};
inline constexpr Unit<-2> iGeV2{1.0e-6, "GeV^-2"}; // 1 GeV^-2 = 1e-6 MeV^-2
inline constexpr Unit<-2> fm2{kFm * kFm, "fm^2"};
inline constexpr Unit<-2> b{kBarn, "b"};
inline constexpr Unit<-2> mb{kBarn * 1.0e-3, "mb"};
inline constexpr Unit<-2> ub{kBarn * 1.0e-6, "ub"};
inline constexpr Unit<-2> nb{kBarn * 1.0e-9, "nb"};
inline constexpr Unit<-2> pb{kBarn * 1.0e-12, "pb"};
inline constexpr Unit<-2> fb{kBarn * 1.0e-15, "fb"};

// ---------------------------------------------------------------------------
// User-defined literals (strongly typed).
// ---------------------------------------------------------------------------
namespace literals {
constexpr Energy operator"" _MeV(long double x) {
    return static_cast<double>(x) * MeV;
}
constexpr Energy operator"" _MeV(unsigned long long x) {
    return static_cast<double>(x) * MeV;
}
constexpr Energy operator"" _GeV(long double x) {
    return static_cast<double>(x) * GeV;
}
constexpr Energy operator"" _GeV(unsigned long long x) {
    return static_cast<double>(x) * GeV;
}
constexpr Length operator"" _fm(long double x) {
    return static_cast<double>(x) * fm;
}
constexpr Length operator"" _fm(unsigned long long x) {
    return static_cast<double>(x) * fm;
}
constexpr Length operator"" _m(long double x) {
    return static_cast<double>(x) * m;
}
constexpr CrossSection operator"" _nb(long double x) {
    return static_cast<double>(x) * nb;
}
constexpr CrossSection operator"" _mb(long double x) {
    return static_cast<double>(x) * mb;
}
} // namespace literals

} // namespace achilles::units

namespace achilles {

// Angle literals. Angles are dimensionless, so they are plain doubles and are
// deliberately outside the Quantity system; a strong Angle type is a separate
// piece of work.
constexpr double operator"" _rad(long double x) {
    return static_cast<double>(x);
}
constexpr double operator"" _deg(long double x) {
    constexpr double ToRads = M_PI / 180;
    return static_cast<double>(x) * ToRads;
}

// The strong literals (_MeV, _GeV, _fm, _m, _mb, _nb) are visible unqualified
// inside namespace achilles, which is where nearly all of the code lives. Unit
// *names* are deliberately not pulled in -- `m`, `b` and `fm` are far too
// common as local variable names -- so spell those units::MeV, units::fm, ...
using namespace units::literals;

} // namespace achilles

#endif
