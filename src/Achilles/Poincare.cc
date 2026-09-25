// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/Poincare.hh"
#include "Achilles/ThreeVector.hh"
#include "Achilles/Utilities.hh"

using namespace achilles;

namespace {
using Vec4 = std::array<double, 4>;

double Dot(const Vec4 &a, const Vec4 &b) noexcept {
    return a[0] * b[0] - (a[1] * b[1] + a[2] * b[2] + a[3] * b[3]);
}
double P2(const Vec4 &a) noexcept {
    return a[1] * a[1] + a[2] * a[2] + a[3] * a[3];
}
double SmallOMCT(const Vec4 &a, const Vec4 &b) noexcept {
    const double mag = std::sqrt(P2(a) * P2(b));
    const double pq = a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
    const double ct = std::min(std::max(pq / mag, -1.), 1.);
    if(ct < 0.) return 1. - ct;
    const Vec4 cross{0., a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3],
                     a[1] * b[2] - a[2] * b[1]};
    const double st = std::sqrt(P2(cross)) / mag;
    const double st2 = st / (2. * std::sqrt((ct + 1) / 2.));
    return 2. * st2 * st2;
}
} // namespace

Poincare::Poincare(const FourVector &v, const double &rsq)
    : m_type(1), m_l(v.ToArray(units::MeV)), m_rsq(rsq > 0. ? rsq : v.M().native()) {}

Poincare::Poincare(const FourVector &v1, const FourVector &v2, int mode)
    : m_type(mode ? 3 : 2), m_l({1., 0., 0., 0.}), m_rsq(1.) {
    if(m_type == 3) {
        m_l = v1.ToArray(units::MeV);
        m_t = v2.ToArray(units::MeV);
        return;
    }
    const double p1 = v1.P().native(), p2 = v2.P().native();
    const Vec4 b{0., v2[1].native() / p2, v2[2].native() / p2, v2[3].native() / p2};
    m_l = Vec4{0., v1[1].native() / p1, v1[2].native() / p1, v1[3].native() / p1};
    const double lb = Dot(m_l, b);
    for(size_t i = 0; i < 4; ++i) m_t[i] = b[i] + m_l[i] * lb;
    double mt(P2(m_t));
    if(mt != 0.)
        for(auto &x : m_t) x *= 1. / std::sqrt(mt);
    size_t l[3]{1, 2, 3};
    double ml[4] = {0., std::abs(m_l[1]), std::abs(m_l[2]), std::abs(m_l[3])};
    if(ml[l[2]] > ml[l[1]]) std::swap<size_t>(l[1], l[2]);
    if(ml[l[1]] > ml[l[0]]) std::swap<size_t>(l[0], l[1]);
    if(ml[l[2]] > ml[l[1]]) std::swap<size_t>(l[1], l[2]);
    double tdp(m_t[l[1]] * m_l[l[1]] + m_t[l[2]] * m_l[l[2]]);
    if(tdp != 0.) m_t[l[0]] = -tdp / m_l[l[0]];
    if(P2(m_t) == 0.) m_t[l[1]] = 1.;
    m_omct = SmallOMCT(m_l, b);
    m_st = -Dot(m_t, b);
}

void Poincare::Boost(FourVector &v) const {
    const auto a = v.ToArray(units::MeV);
    double lv(m_l[1] * a[1] + m_l[2] * a[2] + m_l[3] * a[3]);
    double v0((m_l[0] * a[0] - lv) / m_rsq);
    double c1((a[0] + v0) / (m_rsq + m_l[0]));
    v = FourVector::FromNative({v0, a[1] - c1 * m_l[1], a[2] - c1 * m_l[2], a[3] - c1 * m_l[3]});
}

void Poincare::BoostBack(FourVector &v) const {
    const auto a = v.ToArray(units::MeV);
    double lv(m_l[1] * a[1] + m_l[2] * a[2] + m_l[3] * a[3]);
    double v0((m_l[0] * a[0] + lv) / m_rsq);
    double c1((a[0] + v0) / (m_rsq + m_l[0]));
    v = FourVector::FromNative({v0, a[1] + c1 * m_l[1], a[2] + c1 * m_l[2], a[3] + c1 * m_l[3]});
}

void Poincare::Rotate(FourVector &v) const {
    auto a = v.ToArray(units::MeV);
    double vx(-Dot(m_l, a)), vy(-Dot(m_t, a));
    const double cl = m_omct * vx + m_st * vy;
    const double ct = -m_st * vx + m_omct * vy;
    for(size_t i = 0; i < 4; ++i) a[i] -= cl * m_l[i];
    for(size_t i = 0; i < 4; ++i) a[i] -= ct * m_t[i];
    v = FourVector::FromNative(a);
}

void Poincare::RotateBack(FourVector &v) const {
    auto a = v.ToArray(units::MeV);
    double vx(-Dot(m_l, a)), vy(-Dot(m_t, a));
    const double cl = m_omct * vx - m_st * vy;
    const double ct = m_st * vx + m_omct * vy;
    for(size_t i = 0; i < 4; ++i) a[i] -= cl * m_l[i];
    for(size_t i = 0; i < 4; ++i) a[i] -= ct * m_t[i];
    v = FourVector::FromNative(a);
}

void Poincare::Lambda(FourVector &v) const {
    auto a = v.ToArray(units::MeV);
    const double m2 = Dot(a, a);
    Vec4 lt{};
    for(size_t i = 0; i < 4; ++i) lt[i] = m_l[i] + m_t[i];
    const double c1 = 2.0 * Dot(a, lt) / Dot(lt, lt);
    const double c2 = 2.0 * Dot(a, m_l) / Dot(m_l, m_l);
    for(size_t i = 0; i < 4; ++i) a[i] = a[i] - c1 * lt[i] + c2 * m_t[i];
    a[0] = Sign(a[0]) * std::sqrt(P2(a) + m2);
    v = FourVector::FromNative(a);
}

void Poincare::LambdaBack(FourVector &v) const {
    auto a = v.ToArray(units::MeV);
    const double m2 = Dot(a, a);
    Vec4 lt{};
    for(size_t i = 0; i < 4; ++i) lt[i] = m_l[i] + m_t[i];
    const double c1 = 2.0 * Dot(a, lt) / Dot(lt, lt);
    const double c2 = 2.0 * Dot(a, m_t) / Dot(m_t, m_t);
    for(size_t i = 0; i < 4; ++i) a[i] = a[i] - c1 * lt[i] + c2 * m_l[i];
    a[0] = Sign(a[0]) * std::sqrt(P2(a) + m2);
    v = FourVector::FromNative(a);
}

void Poincare::Invert() {
    if(m_type == 3) {
        std::swap(m_l, m_t);
        return;
    }
    if(m_type == 2) {
        m_st = -m_st;
        return;
    }
    for(size_t i(1); i < 4; ++i) m_l[i] = -m_l[i];
}

FourVector Poincare_Sequence::operator*(const FourVector &p) const {
    FourVector np(p);
    for(const auto &transform : *this) np = transform * np;
    return np;
}

void Poincare_Sequence::Invert() {
    Poincare_Sequence copy(*this);
    reverse_iterator cit(copy.rbegin());
    for(iterator pit(begin()); pit != end(); ++pit, ++cit) {
        cit->Invert();
        *pit = *cit;
    }
}
