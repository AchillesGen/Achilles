// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/CurvedPath.hh"
#include "Achilles/Potential.hh"
#include "Achilles/PotentialFlow.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

using namespace achilles;

// ---------------------------------------------------------------------------------------------
// PairKernel
// ---------------------------------------------------------------------------------------------

double PairKernel::ContactRadius(double sigma) {
    return std::sqrt(sigma / M_PI);
}

double PairKernel::Probability(double b2, double sigma) const {
    switch(m_profile) {
    case Profile::Gaussian:
        return std::exp(-M_PI * b2 / sigma);
    case Profile::Pion:
        return std::exp(-std::sqrt(2 * M_PI / sigma) * std::sqrt(b2));
    case Profile::Cylinder:
        return b2 < sigma / M_PI ? 1 : 0;
    }
    return 0;
}

double PairKernel::Exponent(double b2, double sigma) const {
    const double prob = Probability(b2, sigma);
    if(prob >= 1) return std::numeric_limits<double>::infinity();
    return -std::log1p(-prob);
}

// zeta(1/2 - k), k = 0..15
static constexpr std::array<double, 16> cZetaHalf{
    -1.4603545088095868,   -0.20788622497735457, -0.025485201889833036, 0.008516928777850331,
    0.004441011335479432,  -0.0030916692472158338, -0.0026714580198992244,
    0.0027467679395368687, 0.00326903957260022,   -0.00441603287300489,
    -0.006672172296466641, 0.011146122473942813,  0.02039697871594279,
    -0.04057496748119458,  -0.08717525590621725,  0.2011740493842269};

double PairKernel::LiHalfExp(double t) {
    if(t <= 0) return std::numeric_limits<double>::infinity();
    if(t < 1) {
        // Li_s(e^{-t}) = Gamma(1-s) t^{s-1} + sum_k zeta(s-k) (-t)^k / k!, s = 1/2
        double acc = std::sqrt(M_PI / t);
        double term = 1;
        for(size_t k = 0; k < cZetaHalf.size(); ++k) {
            if(k > 0) term *= -t / static_cast<double>(k);
            acc += cZetaHalf[k] * term;
        }
        return acc;
    }
    double acc = 0;
    for(size_t n = 1; n < 80; ++n) {
        const double dn = static_cast<double>(n);
        const double term = std::exp(-dn * t) / std::sqrt(dn);
        acc += term;
        if(term < 1e-18 * acc) break;
    }
    return acc;
}

double PairKernel::SumK0(double x) {
    if(x <= 0) return std::numeric_limits<double>::infinity();
    if(x < 0.6) {
        // Poisson summation: pi/(2x) + (gamma + ln(x/4pi))/2 + pi sum_l [1/sqrt(x^2+(2 pi l)^2)
        // - 1/(2 pi l)], with the large-l tail in closed form
        constexpr double euler = 0.5772156649015329;
        constexpr size_t lmax = 200;
        double acc = M_PI / (2 * x) + 0.5 * (euler + std::log(x / (4 * M_PI)));
        double sum = 0;
        for(size_t l = 1; l <= lmax; ++l) {
            const double a = 2 * M_PI * static_cast<double>(l);
            sum += 1 / std::sqrt(x * x + a * a) - 1 / a;
        }
        const double dl = static_cast<double>(lmax);
        sum -= x * x / (2 * std::pow(2 * M_PI, 3)) / (2 * dl * dl);
        return acc + M_PI * sum;
    }
    double acc = 0;
    for(size_t n = 1; n < 1000; ++n) {
        const double arg = static_cast<double>(n) * x;
        if(arg > 700) break;
        const double term = std::cyl_bessel_k(0.0, arg);
        acc += term;
        if(term < 1e-18 * acc) break;
    }
    return acc;
}

double PairKernel::Kappa(double r, double sigma) const {
    switch(m_profile) {
    case Profile::Gaussian:
        return LiHalfExp(M_PI * r * r / sigma) / std::sqrt(sigma);
    case Profile::Pion: {
        const double c = std::sqrt(2 * M_PI / sigma);
        return c / M_PI * SumK0(c * r);
    }
    case Profile::Cylinder:
        return r < ContactRadius(sigma) ? std::numeric_limits<double>::infinity() : 0;
    }
    return 0;
}

// ---------------------------------------------------------------------------------------------
// TabulatedPotential
// ---------------------------------------------------------------------------------------------

TabulatedPotential::TabulatedPotential(const Potential &base, double pmax, double dp, double rmax,
                                       double dr)
    : m_np{static_cast<size_t>(std::ceil(pmax / dp)) + 1},
      m_nr{static_cast<size_t>(std::ceil(rmax / dr)) + 1}, m_dp{dp}, m_dr{dr} {
    m_vector.resize(m_np * m_nr);
    m_scalar.resize(m_np * m_nr);
    for(size_t i = 0; i < m_np; ++i) {
        const double p = static_cast<double>(i) * dp;
        for(size_t j = 0; j < m_nr; ++j) {
            const double r = (static_cast<double>(j) + 0.5) * dr;
            const auto vals = base(p, r);
            m_vector[i * m_nr + j] = vals.rvector;
            m_scalar[i * m_nr + j] = vals.rscalar;
        }
    }
}

namespace {
double CatmullRom(double y0, double y1, double y2, double y3, double s) {
    return y1 + 0.5 * s *
                    (y2 - y0 + s * (2 * y0 - 5 * y1 + 4 * y2 - y3 + s * (3 * (y1 - y2) + y3 - y0)));
}
} // namespace

double TabulatedPotential::Interp(const std::vector<double> &table, double p, double r) const {
    r = std::abs(r);
    p = std::clamp(std::abs(p), 0.0, static_cast<double>(m_np - 1) * m_dp);
    // radial coordinate in units of dr, measured from the first node at dr/2
    const double xr = r / m_dr - 0.5;
    const double xp = p / m_dp;
    const auto ir = static_cast<long>(std::floor(xr));
    const auto ip = static_cast<long>(std::floor(xp));
    const double sr = xr - static_cast<double>(ir), sp = xp - static_cast<double>(ip);
    const auto nr = static_cast<long>(m_nr), np = static_cast<long>(m_np);
    auto at = [&](long i, long j) {
        i = std::clamp(i, 0L, np - 1);
        // even extension in r: node -1 mirrors node 0, -2 mirrors 1, ...
        if(j < 0) j = -j - 1;
        j = std::clamp(j, 0L, nr - 1);
        return table[static_cast<size_t>(i) * m_nr + static_cast<size_t>(j)];
    };
    std::array<double, 4> rows{};
    for(long a = 0; a < 4; ++a) {
        const long i = ip - 1 + a;
        rows[static_cast<size_t>(a)] =
            CatmullRom(at(i, ir - 1), at(i, ir), at(i, ir + 1), at(i, ir + 2), sr);
    }
    return CatmullRom(rows[0], rows[1], rows[2], rows[3], sp);
}

namespace {
double CatmullRomDerivative(double y0, double y1, double y2, double y3, double s) {
    return 0.5 * (y2 - y0 + 2 * s * (2 * y0 - 5 * y1 + 4 * y2 - y3) +
                  3 * s * s * (3 * (y1 - y2) + y3 - y0));
}
} // namespace

void TabulatedPotential::Interp3(const std::vector<double> &table, double p, double r,
                                 double *out) const {
    const double sign_r = r < 0 ? -1 : 1;
    r = std::abs(r);
    const double pmaxv = static_cast<double>(m_np - 1) * m_dp;
    const bool clamped = std::abs(p) >= pmaxv;
    p = std::clamp(std::abs(p), 0.0, pmaxv);
    const double xr = r / m_dr - 0.5;
    const double xp = p / m_dp;
    const auto ir = static_cast<long>(std::floor(xr));
    const auto ip = static_cast<long>(std::floor(xp));
    const double sr = xr - static_cast<double>(ir), sp = xp - static_cast<double>(ip);
    const auto nr = static_cast<long>(m_nr), np = static_cast<long>(m_np);
    auto at = [&](long i, long j) {
        i = std::clamp(i, 0L, np - 1);
        if(j < 0) j = -j - 1;
        j = std::clamp(j, 0L, nr - 1);
        return table[static_cast<size_t>(i) * m_nr + static_cast<size_t>(j)];
    };
    std::array<double, 4> rows{}, drows{};
    for(long a = 0; a < 4; ++a) {
        const long i = ip - 1 + a;
        const double y0 = at(i, ir - 1), y1 = at(i, ir), y2 = at(i, ir + 1), y3 = at(i, ir + 2);
        rows[static_cast<size_t>(a)] = CatmullRom(y0, y1, y2, y3, sr);
        drows[static_cast<size_t>(a)] = CatmullRomDerivative(y0, y1, y2, y3, sr) / m_dr;
    }
    out[0] = CatmullRom(rows[0], rows[1], rows[2], rows[3], sp);
    out[1] = clamped ? 0 : CatmullRomDerivative(rows[0], rows[1], rows[2], rows[3], sp) / m_dp;
    out[2] = sign_r * CatmullRom(drows[0], drows[1], drows[2], drows[3], sp);
}

void TabulatedPotential::Evaluate(double p, double r, double *vec, double *scal) const {
    Interp3(m_vector, p, r, vec);
    Interp3(m_scalar, p, r, scal);
}

PotentialVals TabulatedPotential::operator()(double p, double r) const {
    PotentialVals vals{};
    vals.rvector = Interp(m_vector, p, r);
    vals.rscalar = Interp(m_scalar, p, r);
    return vals;
}

// ---------------------------------------------------------------------------------------------
// CurvedPath
// ---------------------------------------------------------------------------------------------

CurvedPath::CurvedPath(double t0, const ThreeVector &q0, const ThreeVector &p0, double mass,
                       std::shared_ptr<Potential> potential, double step, double rexit,
                       bool external)
    : m_t0{t0}, m_mass{mass}, m_step{step}, m_rexit{rexit}, m_external{external},
      m_potential{std::move(potential)} {
    if(m_potential) {
        m_energy = flow::Hamiltonian(q0, p0, *m_potential, m_mass);
        static constexpr double omega = 20;
        m_table = std::dynamic_pointer_cast<TabulatedPotential>(m_potential);
        if(m_table) {
            // same Hamiltonian derivatives as flow::dHdr / flow::dHdp, with the analytic
            // derivatives of the bicubic table in place of five-point stencils
            const auto *table = m_table.get();
            auto dHdr = [mass, table](const ThreeVector &q, const ThreeVector &p,
                                      std::shared_ptr<Potential>) {
                return TableGradient(*table, q, p, mass, true);
            };
            auto dHdp = [mass, table](const ThreeVector &q, const ThreeVector &p,
                                      std::shared_ptr<Potential>) {
                return TableGradient(*table, q, p, mass, false);
            };
            m_integrator = SymplecticIntegrator(q0, p0, m_potential, dHdr, dHdp, omega);
        } else {
            auto dHdr = [mass](const ThreeVector &q, const ThreeVector &p,
                               std::shared_ptr<Potential> pot) {
                return flow::dHdr(q, p, *pot, mass);
            };
            auto dHdp = [mass](const ThreeVector &q, const ThreeVector &p,
                               std::shared_ptr<Potential> pot) {
                return flow::dHdp(q, p, *pot, mass);
            };
            m_integrator = SymplecticIntegrator(q0, p0, m_potential, dHdr, dHdp, omega);
        }
    } else {
        m_energy = std::sqrt(p0.P2() + mass * mass);
    }
    m_nodes.push_back(MakeNode(t0, q0, p0));
}

ThreeVector CurvedPath::TableGradient(const TabulatedPotential &table, const ThreeVector &q,
                                      const ThreeVector &p, double mass, bool radial) {
    const double r = q.P(), pm = p.P();
    double vec[3], scal[3];
    table.Evaluate(pm, r, vec, scal);
    const double meff = mass + scal[0];
    const double energy = std::sqrt(meff * meff + pm * pm);
    if(radial) {
        if(r <= 0) return {0, 0, 0};
        return (meff * scal[2] / energy + vec[2]) * q / r;
    }
    if(pm <= 0) return {0, 0, 0};
    return ((meff * scal[1] + pm) / energy + vec[1]) * p / pm;
}

CurvedPath::Node CurvedPath::MakeNode(double t, const ThreeVector &q, const ThreeVector &p) const {
    Node node{t, q, p, {}, {}};
    if(m_table) {
        node.v = TableGradient(*m_table, q, p, m_mass, false);
        node.f = -TableGradient(*m_table, q, p, m_mass, true);
    } else if(m_potential) {
        node.v = flow::dHdp(q, p, *m_potential, m_mass);
        node.f = -flow::dHdr(q, p, *m_potential, m_mass);
    } else {
        node.v = p / std::sqrt(p.P2() + m_mass * m_mass);
        node.f = {0, 0, 0};
    }
    return node;
}

void CurvedPath::Grow() {
    const size_t k = m_nodes.size();
    const double t = m_t0 + static_cast<double>(k) * m_step;
    if(m_potential) {
        m_integrator.Step<2>(m_step);
        // Tao's two copies (q, p) and (x, y) each carry an O(h^2) ripple that is absent from
        // their average; the averaged state is the node of the dense output (see the note)
        const auto &state = m_integrator.State();
        m_nodes.push_back(MakeNode(t, 0.5 * (state.q + state.x), 0.5 * (state.p + state.y)));
    } else {
        // Exact straight line, evaluated from the start point to avoid accumulated rounding
        const auto &first = m_nodes.front();
        m_nodes.push_back(MakeNode(t, first.q + first.v * (t - m_t0), first.p));
    }
}

const std::vector<CurvedPath::Node> &CurvedPath::Nodes(size_t n) {
    while(m_nodes.size() < n) Grow();
    return m_nodes;
}

size_t CurvedPath::Segment(double t) {
    if(t <= m_t0) {
        Nodes(2);
        return 0;
    }
    // a time within 1e-9 steps below a node belongs to the next segment (rounding of t0 + k h)
    const auto k = static_cast<size_t>(std::floor((t - m_t0) / m_step + 1e-9));
    Nodes(k + 2);
    return k;
}

void CurvedPath::Hermite(size_t k, double t, ThreeVector *q, ThreeVector *v) const {
    const Node &n0 = m_nodes[k];
    const Node &n1 = m_nodes[k + 1];
    const double h = n1.t - n0.t;
    const double s = (t - n0.t) / h;
    const double s2 = s * s, s3 = s2 * s;
    if(q) {
        *q = (2 * s3 - 3 * s2 + 1) * n0.q + (s3 - 2 * s2 + s) * h * n0.v + (-2 * s3 + 3 * s2) * n1.q +
             (s3 - s2) * h * n1.v;
    }
    if(v) {
        *v = ((6 * s2 - 6 * s) * n0.q + (-6 * s2 + 6 * s) * n1.q) / h +
             (3 * s2 - 4 * s + 1) * n0.v + (3 * s2 - 2 * s) * n1.v;
    }
}

ThreeVector CurvedPath::Position(double t) {
    if(Straight()) {
        const auto &first = m_nodes.front();
        return first.q + first.v * (t - m_t0);
    }
    const size_t k = Segment(t);
    ThreeVector q;
    Hermite(k, t, &q, nullptr);
    return q;
}

ThreeVector CurvedPath::Velocity(double t) {
    if(Straight()) return m_nodes.front().v;
    const size_t k = Segment(t);
    ThreeVector v;
    Hermite(k, t, nullptr, &v);
    return v;
}

ThreeVector CurvedPath::Acceleration(double t) {
    if(Straight()) return {0, 0, 0};
    return SegmentAcceleration(Segment(t), t);
}

ThreeVector CurvedPath::SegmentAcceleration(size_t k, double t) {
    if(Straight()) return {0, 0, 0};
    Nodes(k + 2);
    const Node &n0 = m_nodes[k];
    const Node &n1 = m_nodes[k + 1];
    const double h = n1.t - n0.t;
    const double s = (t - n0.t) / h;
    return ((12 * s - 6) * n0.q + (-12 * s + 6) * n1.q) / (h * h) +
           ((6 * s - 4) * n0.v + (6 * s - 2) * n1.v) / h;
}

ThreeVector CurvedPath::Momentum(double t) {
    if(Straight()) return m_nodes.front().p;
    const size_t k = Segment(t);
    const Node &n0 = m_nodes[k];
    const Node &n1 = m_nodes[k + 1];
    const double h = n1.t - n0.t;
    const double s = (t - n0.t) / h;
    const double s2 = s * s, s3 = s2 * s;
    return (2 * s3 - 3 * s2 + 1) * n0.p + (s3 - 2 * s2 + s) * h * n0.f +
           (-2 * s3 + 3 * s2) * n1.p + (s3 - s2) * h * n1.f;
}

double CurvedPath::Hamiltonian(const ThreeVector &q, const ThreeVector &p) const {
    if(!m_potential) return std::sqrt(p.P2() + m_mass * m_mass);
    return flow::Hamiltonian(q, p, *m_potential, m_mass);
}

ThreeVector CurvedPath::ShellMomentum(double t) {
    const auto p = Momentum(t);
    if(Straight()) return p;
    const auto q = Position(t);
    const double k0 = p.P();
    if(k0 <= 0) return p;
    const auto n = p / k0;
    auto f = [&](double k) { return Hamiltonian(q, k * n) - m_energy; };
    double lo = 0, hi = 2 * k0 + 10;
    if(f(lo) > 0 || f(hi) < 0) return p;
    for(size_t it = 0; it < 100 && hi - lo > 1e-10 * hi; ++it) {
        const double mid = 0.5 * (lo + hi);
        (f(mid) > 0 ? hi : lo) = mid;
    }
    return 0.5 * (lo + hi) * n;
}

CurvedPath::Bound CurvedPath::Bounds(double ta, double tb) {
    // the segment holding [ta, tb] (robust against rounding at the end points)
    const size_t k = Segment(0.5 * (ta + tb));
    Nodes(k + 2);
    ThreeVector qa, va, qb, vb;
    if(Straight()) {
        qa = Position(ta);
        qb = Position(tb);
        va = vb = m_nodes.front().v;
        return {qa, qb, 0, va.P(), va.P(), 0, ta, tb};
    }
    Hermite(k, ta, &qa, &va);
    Hermite(k, tb, &qb, &vb);
    const double dt = tb - ta;
    const double amax_k = std::max(SegmentAcceleration(k, ta).P(), SegmentAcceleration(k, tb).P());
    const auto chord = qb - qa;
    // Deviation of a cubic Hermite segment from its chord:
    // q(s) - chord(s) = s(1-s)[(1-s)(dt va - D) - s(dt vb - D)], |.| <= 4/27 (|A| + |C|)
    const double delta = 4.0 / 27.0 * ((dt * va - chord).P() + (dt * vb - chord).P());
    // dq/dt is a quadratic Bezier curve with control points va, 3D/dt - va - vb, vb
    const double vmid = (3 * chord / dt - va - vb).P();
    const double vmax = std::max({va.P(), vb.P(), vmid});
    // d^2q/dt^2 is linear on the segment: its norm is maximal at an end point
    const double amax = amax_k;
    const double vmin = std::max(0.0, std::min(va.P(), vb.P()) - 0.5 * amax * dt);
    return {qa, qb, delta, vmax, vmin, amax, ta, tb};
}

namespace {
double PointSegmentDistance(const achilles::ThreeVector &x, const achilles::ThreeVector &a,
                            const achilles::ThreeVector &b) {
    const auto ab = b - a;
    const double len2 = ab.P2();
    double s = len2 > 0 ? (x - a).Dot(ab) / len2 : 0;
    s = std::clamp(s, 0.0, 1.0);
    return (x - (a + s * ab)).P();
}
} // namespace

const CurvedPath::Bound &CurvedPath::Block(size_t level, size_t index) {
    if(m_blocks.size() <= level) {
        m_blocks.resize(level + 1);
        m_block_set.resize(level + 1);
    }
    auto &blocks = m_blocks[level];
    auto &set = m_block_set[level];
    if(blocks.size() <= index) {
        blocks.resize(index + 1);
        set.resize(index + 1, 0);
    }
    if(set[index]) return blocks[index];
    Bound bound{};
    if(level == 0) {
        const double ta = m_t0 + static_cast<double>(index) * m_step;
        Nodes(index + 2);
        bound = Bounds(ta, m_nodes[index + 1].t);
    } else {
        const Bound left = Block(level - 1, 2 * index);
        const Bound right = Block(level - 1, 2 * index + 1);
        bound.a = left.a;
        bound.b = right.b;
        bound.ta = left.ta;
        bound.tb = right.tb;
        // Every point of a child lies within child.delta of the child's chord, and the child's
        // chord within max(dist(end points)) of this chord (distance to a segment is convex)
        double delta = 0;
        for(const Bound *child : {&left, &right}) {
            const double off = std::max(PointSegmentDistance(child->a, bound.a, bound.b),
                                        PointSegmentDistance(child->b, bound.a, bound.b));
            delta = std::max(delta, child->delta + off);
        }
        bound.delta = delta;
        bound.vmax = std::max(left.vmax, right.vmax);
        bound.vmin = std::min(left.vmin, right.vmin);
        bound.amax = std::max(left.amax, right.amax);
    }
    // the vectors may have been reallocated by the recursion
    m_blocks[level][index] = bound;
    m_block_set[level][index] = 1;
    return m_blocks[level][index];
}

size_t CurvedPath::ExitSegment() {
    const double texit = ExitTime();
    if(!std::isfinite(texit)) return std::numeric_limits<size_t>::max();
    return Segment(texit);
}

double CurvedPath::MaxSpeed(double tmax) {
    if(Straight()) return m_nodes.front().v.P();
    double vmax = 0;
    for(size_t k = 0;; ++k) {
        const double ta = m_t0 + static_cast<double>(k) * m_step;
        if(ta >= tmax) break;
        vmax = std::max(vmax, Bounds(ta, ta + m_step).vmax);
    }
    return vmax;
}

bool CurvedPath::Trapped() {
    // Central field: E and L = |q x p| are conserved, and the radial motion is confined to the
    // connected region around r(t0) where k(r)^2 r^2 >= L^2, k(r) solving H(k, r) = E.
    // The particle can leave only if that region reaches the exit radius.
    if(!m_potential) return false;
    const auto &first = m_nodes.front();
    const double r0 = first.q.P();
    if(r0 >= m_rexit) return false;
    const double ell2 = first.q.Cross(first.p).P2();
    auto allowed = [&](double r) {
        // k(r)^2 r^2 >= L^2  <=>  H(L/r, r) <= E, since H increases with |p|
        return flow::Hamiltonian({r, 0, 0}, {std::sqrt(ell2) / r, 0, 0}, *m_potential, m_mass) <=
               m_energy;
    };
    constexpr double dr = 0.01;
    for(double r = r0 + dr; r < m_rexit + dr; r += dr) {
        if(!allowed(std::min(r, m_rexit))) return true;
    }
    return false;
}

double CurvedPath::ExitTime() {
    if(m_exit >= 0) return m_exit;
    const double r2 = m_rexit * m_rexit;
    auto outside = [&](const ThreeVector &q) {
        return q.P2() > r2 && (!m_external || q.Pz() >= m_rexit);
    };
    if(!m_external && outside(m_nodes.front().q)) return m_exit = m_t0;
    if(Straight()) {
        // Same root as Cascade::ExitTime, without its margin
        const auto &first = m_nodes.front();
        const auto pos = first.q, vel = first.v;
        const double a = vel.P2();
        if(a <= 0) return m_exit = std::numeric_limits<double>::infinity();
        const double b = pos.Dot(vel), c = pos.P2() - r2;
        double time = std::max(0.0, (-b + std::sqrt(std::max(0.0, b * b - a * c))) / a);
        if(m_external) {
            if(vel.Pz() <= 0) return m_exit = std::numeric_limits<double>::infinity();
            time = std::max(time, (m_rexit - pos.Pz()) / vel.Pz());
        }
        return m_exit = m_t0 + time;
    }
    static constexpr size_t max_segments = 50000;
    for(size_t k = 0; k < max_segments; ++k) {
        Nodes(k + 2);
        if(!outside(m_nodes[k + 1].q)) continue;
        double lo = m_nodes[k].t, hi = m_nodes[k + 1].t;
        for(size_t it = 0; it < 80; ++it) {
            const double mid = 0.5 * (lo + hi);
            ThreeVector q;
            Hermite(k, mid, &q, nullptr);
            (outside(q) ? hi : lo) = mid;
        }
        return m_exit = hi;
    }
    return m_exit = std::numeric_limits<double>::infinity();
}
