// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/FourVector.hh"
#include "Achilles/ThreeVector.hh"

extern "C" {

// The Fortran side works in MeV for momenta and fm for positions; every value
// crossing this boundary is converted explicitly.
using Vec4 = achilles::FourVector;
using Vec3 = achilles::ThreeVector;
using Boost3 = achilles::ThreeBoost;
namespace units = achilles::units;

// FourVector Constructor and Destructor
Vec4 *CreateFourVector(const double E, const double x, const double y, const double z) {
    return new Vec4(Vec4::FromNative({E, x, y, z}));
}
Vec4 *CopyFourVector(const Vec4 *self) {
    return new Vec4(*self);
}
void DeleteFourVector(Vec4 *self) {
    delete self;
}

// Functions with FourVectors
Vec4 *Boost(const Vec4 *self, const Vec3 *beta) {
    // The Fortran side carries a velocity in a plain triple; it is dimensionless.
    const auto b = beta->Native();
    return new Vec4(
        self->Boost(Boost3{units::Dimensionless{b[0].native()}, units::Dimensionless{b[1].native()},
                           units::Dimensionless{b[2].native()}}));
}
Vec3 *BoostVector(const Vec4 *self) {
    const auto b = self->BoostVector().Native();
    return new Vec3(Vec3::FromNative({b[0].native(), b[1].native(), b[2].native()}));
}
double Dot4(const Vec4 *self, const Vec4 *other) {
    return self->Dot(*other).native();
}
Vec4 *Add4(const Vec4 *self, const Vec4 *other) {
    return new Vec4((*self) + (*other));
}
Vec4 *Sub4(const Vec4 *self, const Vec4 *other) {
    return new Vec4((*self) - (*other));
}
Vec4 *Scale4(const Vec4 *self, const double scale) {
    return new Vec4(scale * (*self));
}
double Get4(const Vec4 *self, const int idx) {
    return self->operator[](static_cast<size_t>(idx)).native();
}
void Print4(const Vec4 *self) {
    fmt::print("{}\n", *self);
}

// ThreeVector Constructor and Destructor
Vec3 *CreateThreeVector(const double x, const double y, const double z) {
    return new Vec3(Vec3::FromNative({x, y, z}));
}
Vec3 *New3(Vec3 *self) {
    return new Vec3(*self);
}

void DeleteThreeVector(Vec3 *self) {
    delete self;
}

// Functions with ThreeVectors
double Dot3(const Vec3 *self, const Vec3 *other) {
    return self->Dot(*other).native();
}
Vec3 *Add3(const Vec3 *self, const Vec3 *other) {
    return new Vec3((*self) + (*other));
}
Vec3 *Sub3(const Vec3 *self, const Vec3 *other) {
    return new Vec3((*self) - (*other));
}
Vec3 *Scale3(const Vec3 *self, const double scale) {
    return new Vec3(scale * (*self));
}
double Get3(const Vec3 *self, const int idx) {
    return self->operator[](static_cast<size_t>(idx)).native();
}
void Print3(const Vec3 *self) {
    fmt::print("{}\n", *self);
}
}
