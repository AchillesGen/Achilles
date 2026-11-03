// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

// Curved-path, event-driven cascade algorithms on the frozen background of an event:
//
//  CurvedVeto (option a): one Bernoulli trial P(b^2, sigma) per local minimum of the pair
//      distance D(t) on the dense path, i.e. Algorithm: Veto generalised to curved paths.
//      The tube bound P(d_low^2, sigma) >= P on each integration step is the overestimate; a
//      single uniform decides both the prescreen and the acceptance, so a minimum is located
//      only when the prescreen passes.
//
//  Hazard (option b): each (propagating, background) pair is a channel with hazard
//      Gamma(t) = |dq/dt| kappa(|q(t) - r_j|; sigma), kappa the inverse Abel transform of
//      X(b) = -ln(1 - P(b^2, sigma)). Candidates come from a piecewise-constant majorant
//      vmax kappa(d_low) on (sub)steps of the dense path (inverted exactly), are accepted
//      with Gamma / Gamma_hat at the candidate, and a rejected candidate continues the same
//      channel from that time. All channels compete in one time-ordered queue; the earliest
//      accepted candidate is tried first. Pauli blocking thins the accepted candidates.
//
// The flow between collisions is Hamilton's equations in the static, real nuclear mean field
// (nucleons only, PotentialProp: True) integrated with SymplecticIntegrator; without a
// potential the paths are exact straight lines and CurvedVeto reproduces Algorithm: Veto.

#include "Achilles/Cascade.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Event.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Potential.hh"
#include "Achilles/PotentialFlow.hh"
#include "Achilles/Random.hh"

#include "spdlog/spdlog.h"

#include <algorithm>
#include <cmath>
#include <array>
#include <functional>
#include <limits>

using namespace achilles;

namespace {

constexpr double cInf = std::numeric_limits<double>::infinity();

/// Distance from x to the segment [a, b]
double SegmentDistance(const ThreeVector &x, const ThreeVector &a, const ThreeVector &b) {
    const auto ab = b - a;
    const double len2 = ab.P2();
    double s = len2 > 0 ? (x - a).Dot(ab) / len2 : 0;
    s = std::clamp(s, 0.0, 1.0);
    return (x - (a + s * ab)).P();
}

double Exponential() {
    return -std::log(Random::Instance().Uniform(0.0, 1.0));
}

/// Zero-truncated Poisson variate with mean parameter x
size_t TruncatedPoisson(double x) {
    if(!(x < 700)) return 1000;
    // inversion on k >= 1
    const double norm = -std::expm1(-x);
    const double u = Random::Instance().Uniform(0.0, 1.0) * norm;
    double p = std::exp(-x) * x;
    double cdf = p;
    size_t k = 1;
    while(cdf < u && k < 1000) {
        ++k;
        p *= x / static_cast<double>(k);
        cdf += p;
    }
    return k;
}

} // namespace

Cascade::CurvedOptions::Pauli Cascade::PauliMode() const {
    if(m_curved.pauli != CurvedOptions::Pauli::Default) return m_curved.pauli;
    return m_algorithm == Algorithm::Hazard ? CurvedOptions::Pauli::Thinning
                                            : CurvedOptions::Pauli::SingleAttempt;
}

std::shared_ptr<Potential> Cascade::CurvedPotential() const {
    auto base = m_nucleus->GetPotential();
    if(!base || !m_curved.tabulate) return base;
    if(m_table_source != base.get()) {
        m_table = std::make_shared<TabulatedPotential>(*base);
        m_table_source = base.get();
    }
    return m_table;
}

void Cascade::RecordFirst(const Particles &particles, size_t idx1, size_t idx2) {
    if(m_diag.first_set || m_first_kicked.count(idx1) == 0) return;
    m_diag.first_set = true;
    m_diag.first_time = currentTime;
    m_diag.first_partner = idx2;
    m_diag.first_position = particles[idx1].Position();
    m_diag.first_separation = (particles[idx1].Position() - particles[idx2].Position()).P();
    const double v = m_first_v0.P();
    if(v > 0) {
        const auto d = particles[idx2].Position() - m_first_q0;
        const auto vhat = m_first_v0 / v;
        m_diag.first_z = d.Dot(vhat);
        m_diag.first_b = (d - m_diag.first_z * vhat).P();
    }
}

std::optional<Cascade::Track> Cascade::MakeTrack(Event &event, size_t idx, bool with_potential,
                                                 double start) const {
    const Particle &part = event.Hadrons()[idx];
    const bool pot = with_potential && part.Info().IsNucleon();
    const double mass = part.Info().Mass();
    ThreeVector p0 = part.Momentum().Vec3();
    const bool external = part.Status() == ParticleStatus::external_test;
    std::shared_ptr<Potential> potential = pot ? CurvedPotential() : nullptr;
    if(pot && external && p0.P() > 0) {
        // A beam particle is defined by its asymptotic momentum: put it on the shell
        // H(q0, p) = sqrt(p_inf^2 + m^2) at its starting point
        const double energy = std::sqrt(p0.P2() + mass * mass);
        const auto n = p0 / p0.P();
        double lo = 0, hi = 2 * p0.P() + 100;
        for(size_t it = 0; it < 200; ++it) {
            const double mid = 0.5 * (lo + hi);
            (flow::Hamiltonian(part.Position(), mid * n, *potential, mass) > energy ? hi : lo) =
                mid;
        }
        p0 = 0.5 * (lo + hi) * n;
    }
    Track track{CurvedPath(start, part.Position(), p0, mass, potential, m_curved.path_step,
                           m_nucleus->Radius(), external),
                start, pot};
    if(part.InFormationZone() && !part.Info().IsPion()) track.t_formed += part.FormationZone();
    if(pot && (track.path.Energy() < mass || track.path.Trapped())) return std::nullopt;
    return track;
}

bool Cascade::StartTrack(Event &event, size_t idx, bool with_potential) {
    Particles &particles = event.Hadrons();
    auto track = MakeTrack(event, idx, with_potential, currentTime);
    if(!track) {
        // Bound in the mean field: captured where it was created
        particles[idx].Status() = ParticleStatus::captured;
        ++m_diag.captured;
        return false;
    }
    if(track->potential) {
        const auto p = track->path.Momentum(currentTime);
        const double mass = track->path.Mass();
        particles[idx].SetMomentum({std::sqrt(p.P2() + mass * mass), p.Px(), p.Py(), p.Pz()});
    }
    const double texit = track->path.ExitTime();
    if(!std::isfinite(texit)) {
        // Not trapped by the barrier test but no exit within the time cap of the path
        particles[idx].Status() = ParticleStatus::captured;
        ++m_diag.captured;
        return false;
    }
    m_tracks[idx] = std::move(*track);
    if(std::isfinite(texit)) m_time_steps.push({texit, {idx, cExitTrial}});
    if(!particles[idx].Info().IsStable()) ScheduleDecay(particles[idx], idx);
    return true;
}

double Cascade::PairSigma(Event &event, size_t idx, size_t hit) {
    const auto key = std::make_pair(idx, hit);
    auto it = m_pair_sigma.find(key);
    if(it != m_pair_sigma.end()) return it->second;
    // Asymptotic momentum of the propagating particle: |p_inf| = sqrt(E^2 - m^2), constant
    // between collisions, along its current direction
    Particle &part = event.Hadrons()[idx];
    const auto saved = part.Momentum();
    const auto &track = m_tracks.at(idx);
    if(track.potential) {
        const double mass = track.path.Mass();
        const double energy = track.path.Energy();
        const auto p = saved.Vec3();
        const double pinf = std::sqrt(std::max(0.0, energy * energy - mass * mass));
        const auto pvec = p.P() > 0 ? pinf * p / p.P() : p;
        part.SetMomentum({energy, pvec.Px(), pvec.Py(), pvec.Pz()});
    }
    double sigma = GetXSec(event, idx, hit) / 10;
    part.SetMomentum(saved);
    if(!(sigma > 0)) sigma = 0;
    m_pair_sigma[key] = sigma;
    return sigma;
}

void Cascade::ScheduleCurvedPairs(Event &event, size_t idx) {
    const Particles &particles = event.Hadrons();
    for(size_t j = 0; j < particles.size(); ++j) {
        if(particles[j].Status() != ParticleStatus::background) continue;
        SchedulePair(event, idx, j, currentTime);
    }
}

void Cascade::SchedulePair(Event &event, size_t idx, size_t hit, double start) {
    const double sigma = PairSigma(event, idx, hit);
    if(!(sigma > 0)) return;
    const ThreeVector x = event.Hadrons()[hit].Position();
    m_event = &event;
    m_chan_idx = idx;
    m_chan_hit = hit;
    m_local_sigma = m_curved.sigma == CurvedOptions::Sigma::Local;
    std::optional<double> time;
    if(m_algorithm == Algorithm::Hazard) {
        auto &track = m_tracks.at(idx);
        const double from = std::max(start, track.t_formed);
        time = m_kernel.HardSphere() ? NextContactTime(idx, x, sigma, from)
                                     : NextHazardTime(idx, x, sigma, from);
    } else {
        time = NextMinimumTrial(idx, hit, x, sigma, start);
    }
    m_local_sigma = false;
    if(time) m_time_steps.push({*time, {idx, hit}});
}

double Cascade::LocalSigma(double t) {
    // Pair cross section with the local (dense-output) momentum of the propagating particle
    Particle &part = m_event->Hadrons()[m_chan_idx];
    const auto saved = part.Momentum();
    const auto p = m_tracks.at(m_chan_idx).path.Momentum(t);
    const double mass = part.Info().Mass();
    part.SetMomentum({std::sqrt(p.P2() + mass * mass), p.Px(), p.Py(), p.Pz()});
    const double sigma = GetXSec(*m_event, m_chan_idx, m_chan_hit) / 10;
    part.SetMomentum(saved);
    return sigma > 0 ? sigma : 0;
}

double Cascade::KappaBound(double r, double sigma_lo, double sigma_hi) const {
    // sup of kappa(r; sigma) for sigma in [sigma_lo, sigma_hi]
    switch(m_kernel.Type()) {
    case PairKernel::Profile::Gaussian:
        return PairKernel::LiHalfExp(M_PI * r * r / sigma_hi) / std::sqrt(sigma_lo);
    case PairKernel::Profile::Pion: {
        const double chi = std::sqrt(2 * M_PI / sigma_lo), clo = std::sqrt(2 * M_PI / sigma_hi);
        return chi / M_PI * PairKernel::SumK0(clo * r);
    }
    case PairKernel::Profile::Cylinder:
        return m_kernel.Kappa(r, sigma_hi);
    }
    return cInf;
}

double Cascade::PairHazard(Track &track, const ThreeVector &x, double sigma, double t) {
    if(m_local_sigma) sigma = LocalSigma(t);
    if(!(sigma > 0)) return 0;
    return track.path.Velocity(t).P() * m_kernel.Kappa((track.path.Position(t) - x).P(), sigma);
}

double Cascade::HazardMajorant(Track &track, const ThreeVector &x, double sigma, double ta,
                               double *tb) {
    // Piecewise-constant majorant vmax kappa(d_low) on the next interval [ta, tb] of the dense
    // path, d_low a lower bound on the pair distance over the interval. Aligned with a node,
    // the interval is the largest block of 2^L segments whose majorant exponent is at most
    // max_step_exponent (blocks are shared by all pairs of the path); otherwise it is the rest
    // of the current segment, halved until the exponent is small enough.
    auto &path = track.path;
    const double hmax = m_curved.max_step_exponent;
    auto majorant = [&](const CurvedPath::Bound &bound) {
        const double dlow = SegmentDistance(x, bound.a, bound.b) - bound.delta;
        if(!(dlow > 0)) return cInf;
        if(!m_local_sigma) return bound.vmax * m_kernel.Kappa(dlow, sigma);
        // local cross section: range from three samples widened by 15% (a violation, if any,
        // is caught at the candidate and handled by the weighted veto)
        const double s1 = LocalSigma(bound.ta), s2 = LocalSigma(bound.tb),
                     s3 = LocalSigma(0.5 * (bound.ta + bound.tb));
        const double lo = 0.85 * std::min({s1, s2, s3}), hi = 1.15 * std::max({s1, s2, s3});
        if(!(hi > 0)) return 0.0;
        return bound.vmax * KappaBound(dlow, std::max(lo, 1e-6), hi);
    };
    const size_t k = path.Segment(ta);
    const double tk = path.SegmentStart(k);
    if(std::abs(ta - tk) <= 1e-12 * std::max(1.0, std::abs(tk))) {
        const size_t kexit = path.ExitSegment();
        for(size_t level = CurvedPath::cMaxLevel + 1; level-- > 1;) {
            const size_t width = size_t{1} << level;
            if(k % width != 0 || k + width > kexit + 1) continue;
            const auto &bound = path.Block(level, k / width);
            const double gmax = majorant(bound);
            if(gmax * (bound.tb - bound.ta) <= hmax) {
                *tb = bound.tb;
                return gmax;
            }
        }
    }
    double tend = path.SegmentEnd(k);
    double gmax = cInf;
    for(size_t depth = 0;; ++depth) {
        gmax = majorant(path.Bounds(ta, tend));
        if(gmax * (tend - ta) <= hmax || depth >= 60 || tend - ta < 1e-12) break;
        tend = ta + 0.5 * (tend - ta);
        ++m_diag.subdivisions;
    }
    *tb = tend;
    return gmax;
}

std::optional<double> Cascade::NextHazardTime(size_t idx, const ThreeVector &x, double sigma,
                                              double start) {
    auto &track = m_tracks.at(idx);
    const double texit = track.path.ExitTime();
    if(start >= texit) return std::nullopt;

    // Acceptance of a candidate at time tc against the majorant value gmax. Returns true to
    // accept. A violated majorant (Gamma > Gamma_hat) is never clipped: Lonnblad's weighted
    // veto with acceptance 1/2 keeps the expectation exact at the price of event weights.
    auto accept = [&](double tc, double gmax) {
        ++m_diag.candidates;
        const double ratio = PairHazard(track, x, sigma, tc) / gmax;
        m_diag.max_ratio = std::max(m_diag.max_ratio, ratio);
        if(ratio > 1) {
            ++m_diag.violations;
            constexpr double eps = 0.5;
            if(Random::Instance().Uniform(0.0, 1.0) < eps) {
                m_diag.weight *= ratio / eps;
                ++m_diag.accepted;
                return true;
            }
            m_diag.weight *= (1 - ratio) / (1 - eps);
            return false;
        }
        if(Random::Instance().Uniform(0.0, 1.0) < ratio) {
            ++m_diag.accepted;
            return true;
        }
        return false;
    };

    if(m_curved.majorant == CurvedOptions::Majorant::Constant) {
        double gconst = 0;
        for(double ta = start; ta < texit;) {
            double tb{};
            gconst = std::max(gconst, HazardMajorant(track, x, sigma, ta, &tb));
            ta = tb;
        }
        if(!(gconst > 0)) return std::nullopt;
        if(!std::isfinite(gconst)) {
            ++m_diag.violations;
            return start;
        }
        for(double t = start;;) {
            t += Exponential() / gconst;
            if(t >= texit) return std::nullopt;
            if(accept(t, gconst)) return t;
        }
    }

    const double inflate = m_curved.majorant == CurvedOptions::Majorant::Inflated ? 2 : 1;
    double t = start;
    double budget = Exponential();
    while(t < texit) {
        double tb{};
        const double gmax = inflate * HazardMajorant(track, x, sigma, t, &tb);
        tb = std::min(tb, texit);
        if(!std::isfinite(gmax)) {
            // The dense path runs through the spectator (to 2^-60 of a step): the hazard
            // integral diverges there and the pair collides with certainty
            ++m_diag.violations;
            ++m_diag.accepted;
            return t;
        }
        const double hmax = gmax * (tb - t);
        if(hmax < budget) {
            budget -= hmax;
            t = tb;
            continue;
        }
        const double tc = t + budget / gmax;
        if(accept(tc, gmax)) return tc;
        // Rejected: continue the same channel from the candidate (memorylessness)
        t = tc;
        budget = Exponential();
    }
    return std::nullopt;
}

std::optional<double> Cascade::NextContactTime(size_t idx, const ThreeVector &x, double sigma,
                                               double start) {
    auto &track = m_tracks.at(idx);
    auto &path = track.path;
    const double texit = path.ExitTime();
    const double radius = PairKernel::ContactRadius(sigma);
    if(start >= texit) return std::nullopt;
    if((path.Position(start) - x).P() <= radius) {
        ++m_diag.contacts;
        return start;
    }
    for(double t = start; t < texit;) {
        const size_t k = path.Segment(t);
        double tend = path.SegmentEnd(k);
        if(tend <= t) tend = path.SegmentEnd(k + 1);
        tend = std::min(tend, texit);
        const auto bound = path.Bounds(t, tend);
        if(SegmentDistance(x, bound.a, bound.b) - bound.delta > radius) {
            t = tend;
            continue;
        }
        constexpr size_t nsample = 32;
        double prev = t;
        for(size_t i = 1; i <= nsample; ++i) {
            const double ti = t + (tend - t) * static_cast<double>(i) / nsample;
            if((path.Position(ti) - x).P() <= radius) {
                double lo = prev, hi = ti;
                for(size_t it = 0; it < 60; ++it) {
                    const double mid = 0.5 * (lo + hi);
                    ((path.Position(mid) - x).P() <= radius ? hi : lo) = mid;
                }
                ++m_diag.contacts;
                return hi;
            }
            prev = ti;
        }
        t = tend;
    }
    return std::nullopt;
}

std::vector<Cascade::Bracket> Cascade::CurvedMinima(Track &track, const ThreeVector &x,
                                                    double start, double stop) {
    // Brackets [ta, tb] with g(ta) < 0 <= g(tb), g = (q - x).dq/dt, on pieces where g is
    // monotone (g' >= vmin^2 - D_max amax > 0), so that each bracket holds exactly one local
    // minimum of the pair distance. Blocks of segments are tested before their children.
    std::vector<Bracket> brackets;
    auto &path = track.path;
    if(path.Straight()) {
        const auto v = path.Velocity(start);
        const double v2 = v.P2();
        if(v2 <= 0) return brackets;
        const double tc = start + (x - path.Position(start)).Dot(v) / v2;
        if(tc > start && tc < stop) brackets.push_back({tc, tc, (path.Position(tc) - x).P()});
        return brackets;
    }
    auto g = [&](double t) { return (path.Position(t) - x).Dot(path.Velocity(t)); };
    auto monotone = [&](const CurvedPath::Bound &bound) {
        const double dmax = std::max((bound.a - x).P(), (bound.b - x).P()) + bound.delta;
        return bound.vmin > 0 && bound.vmin * bound.vmin > dmax * bound.amax;
    };
    // pieces in time order: (level, index) blocks or explicit [ta, tb] intervals
    struct Piece {
        double ta, tb;
        size_t level, index;
        bool block;
    };
    std::vector<Piece> work; // used as a stack, last = earliest
    std::vector<Piece> pieces;
    {
        size_t k = path.Segment(start);
        double t = start;
        if(path.SegmentStart(k) < t) {
            pieces.push_back({t, std::min(path.SegmentEnd(k), stop), 0, k, false});
            t = path.SegmentEnd(k);
            ++k;
        }
        const size_t kstop = path.Segment(stop);
        while(t < stop && k <= kstop) {
            size_t level = CurvedPath::cMaxLevel;
            while(level > 0 && (k % (size_t{1} << level) != 0 || k + (size_t{1} << level) > kstop))
                --level;
            const size_t width = size_t{1} << level;
            if(level == 0 && path.SegmentEnd(k) > stop) {
                pieces.push_back({t, stop, 0, k, false});
                break;
            }
            pieces.push_back({t, path.SegmentEnd(k + width - 1), level, k / width, true});
            t = path.SegmentEnd(k + width - 1);
            k += width;
        }
    }
    work.assign(pieces.rbegin(), pieces.rend());
    while(!work.empty()) {
        const Piece piece = work.back();
        work.pop_back();
        if(piece.tb <= piece.ta) continue;
        const auto bound = piece.block ? path.Block(piece.level, piece.index)
                                       : path.Bounds(piece.ta, piece.tb);
        // Exclusion: with vbar the mean velocity of the piece, |v - vbar| <= amax dt, and every
        // point within delta of the chord, g = (q - x).v has a fixed sign when the projection
        // of x on the chord lies beyond either end by more than the margin
        {
            const double dt = bound.tb - bound.ta;
            const auto chord = bound.b - bound.a;
            const double len = chord.P();
            if(len > 0 && dt > 0) {
                const auto u = chord / len;
                const double vbar = len / dt;
                const double dmax = std::max((bound.a - x).P(), (bound.b - x).P()) + bound.delta;
                const double margin = bound.delta + dmax * bound.amax * dt / vbar;
                const double sx = (x - bound.a).Dot(u);
                if(sx < -margin || sx > len + margin) continue;
            }
        }
        if(!monotone(bound)) {
            if(piece.block && piece.level > 0) {
                work.push_back({0, 0, piece.level - 1, 2 * piece.index + 1, true});
                work.push_back({0, 0, piece.level - 1, 2 * piece.index, true});
                work.back().ta = path.Block(piece.level - 1, 2 * piece.index).ta;
                work.back().tb = path.Block(piece.level - 1, 2 * piece.index).tb;
                work[work.size() - 2].ta = path.Block(piece.level - 1, 2 * piece.index + 1).ta;
                work[work.size() - 2].tb = path.Block(piece.level - 1, 2 * piece.index + 1).tb;
                continue;
            }
            if(piece.tb - piece.ta > path.Step() / 64) {
                const double mid = 0.5 * (piece.ta + piece.tb);
                work.push_back({mid, piece.tb, 0, 0, false});
                work.push_back({piece.ta, mid, 0, 0, false});
                continue;
            }
        }
        if(!monotone(bound)) {
            // below 1/64 of a step without a monotonicity proof: sample g at 9 points (two
            // minima of the pair distance within such a piece need a loop of ~h/64 radius)
            ++m_diag.unresolved;
            double prev = piece.ta;
            for(size_t i = 1; i <= 8; ++i) {
                const double ti = piece.ta + (piece.tb - piece.ta) * static_cast<double>(i) / 8;
                if(g(prev) < 0 && g(ti) >= 0) {
                    const auto sub = path.Bounds(prev, ti);
                    brackets.push_back(
                        {prev, ti, std::max(0.0, SegmentDistance(x, sub.a, sub.b) - sub.delta)});
                }
                prev = ti;
            }
            continue;
        }
        if(g(piece.ta) < 0 && g(piece.tb) >= 0) {
            const double dlow =
                std::max(0.0, SegmentDistance(x, bound.a, bound.b) - bound.delta);
            brackets.push_back({piece.ta, piece.tb, dlow});
        }
    }
    return brackets;
}

std::optional<double> Cascade::NextMinimumTrial(size_t idx, size_t hit, const ThreeVector &x,
                                                double sigma, double start) {
    auto &track = m_tracks.at(idx);
    auto &path = track.path;
    const double texit = path.ExitTime();
    if(start >= texit) return std::nullopt;
    for(const auto &[ta, tb, dlow] : CurvedMinima(track, x, start, texit)) {
        ++m_diag.minima;
        // Tube-bound overestimate on the bracket: P_hat = P(d_low^2, sigma) >= P(b^2, sigma)
        double sigma_hat = sigma;
        if(m_local_sigma) {
            sigma_hat = 1.15 * std::max({LocalSigma(ta), LocalSigma(tb), LocalSigma(0.5 * (ta + tb))});
        }
        const double phat = m_kernel.Probability(dlow * dlow, sigma_hat);
        const double u = Random::Instance().Uniform(0.0, 1.0);
        ++m_diag.candidates;
        if(u >= phat) {
            ++m_diag.prescreened;
            continue;
        }
        // Locate the minimum (bisection on g) only for candidates that pass the prescreen
        double tmin = ta;
        if(tb > ta) {
            double lo = ta, hi = tb;
            for(size_t it = 0; it < 100 && hi - lo > 1e-14 * std::max(1.0, std::abs(hi)); ++it) {
                const double mid = 0.5 * (lo + hi);
                ((path.Position(mid) - x).Dot(path.Velocity(mid)) < 0 ? lo : hi) = mid;
            }
            tmin = 0.5 * (lo + hi);
        }
        // No trial inside the formation zone (pions exempt via t_formed), as in Veto
        if(tmin < track.t_formed) continue;
        const double b2 = (path.Position(tmin) - x).P2();
        const double sigma_min = m_local_sigma ? LocalSigma(tmin) : sigma;
        const double prob = m_kernel.Probability(b2, sigma_min);
        if(prob > phat * (1 + 1e-9)) ++m_diag.violations;
        if(!(u < prob)) continue;
        ++m_diag.accepted;
        if(PauliMode() == CurvedOptions::Pauli::Thinning) {
            // Thinned atom: Poisson(X) attempts, at least one given that it fired
            const size_t attempts = TruncatedPoisson(m_kernel.Exponent(b2, sigma_min));
            m_atom_attempts[{idx, hit}] = attempts - 1;
        }
        return tmin;
    }
    return std::nullopt;
}

void Cascade::SyncParticles(Particles &particles, double time) {
    const double dt = time - currentTime;
    for(size_t i = 0; i < particles.size(); ++i) {
        auto &part = particles[i];
        if(!part.IsPropagating()) continue;
        auto it = m_tracks.find(i);
        if(it == m_tracks.end()) {
            part.Propagate(dt);
            continue;
        }
        auto &path = it->second.path;
        part.SetPosition(path.Position(time));
        if(it->second.potential) {
            const auto p = path.Momentum(time);
            const double mass = path.Mass();
            part.SetMomentum({std::sqrt(p.P2() + mass * mass), p.Px(), p.Py(), p.Pz()});
        }
        if(part.InFormationZone()) part.UpdateFormationZone(dt);
    }
}

void Cascade::FinishTrack(Particles &particles, size_t idx) {
    auto &part = particles[idx];
    auto &track = m_tracks.at(idx);
    if(!track.potential) {
        // Straight line: same classification as Escaped in Algorithm: Veto
        const double speed = track.path.Velocity(currentTime).P();
        const double margin = speed > 0 ? 1e-6 / speed : 0;
        part.SetPosition(track.path.Position(currentTime + margin));
        Escaped(particles);
        if(part.IsPropagating()) {
            // rounding at the sphere: classify with the same rule
            constexpr double potential = 10.0;
            const double energy = part.Momentum().E() - Constant::mN - potential;
            part.Status() = energy > 0 || !part.Info().IsNucleon() ? ParticleStatus::final_state
                                                                  : ParticleStatus::captured;
        }
        return;
    }
    // Mean field: the energy H is conserved; the particle leaves with the asymptotic
    // momentum sqrt(H^2 - m^2) along its direction at the exit radius (potential counted once)
    const double mass = track.path.Mass();
    const double energy = track.path.Energy();
    const auto p = track.path.Momentum(currentTime);
    part.SetPosition(track.path.Position(currentTime));
    if(energy <= mass || p.P() <= 0) {
        part.Status() = ParticleStatus::captured;
        return;
    }
    const double pinf = std::sqrt(energy * energy - mass * mass);
    const auto pvec = pinf * p / p.P();
    part.SetMomentum({energy, pvec.Px(), pvec.Py(), pvec.Pz()});
    part.Status() = ParticleStatus::final_state;
}

void Cascade::EvolveCurved(Event &event, Nucleus *nucleus) {
    m_nucleus = nucleus;
    currentTime = 0;
    m_time_steps = {};
    m_tracks.clear();
    m_pair_sigma.clear();
    m_atom_attempts.clear();
    Particles &particles = event.Hadrons();
    const bool with_potential = m_potential_prop && m_nucleus->GetPotential() != nullptr;

    const size_t ninitial = particles.size();
    for(size_t idx = 0; idx < ninitial; ++idx) {
        if(!particles[idx].IsPropagating()) continue;
        if(StartTrack(event, idx, with_potential)) ScheduleCurvedPairs(event, idx);
    }

    while(!m_time_steps.empty()) {
        const auto entry = m_time_steps.top();
        m_time_steps.pop();
        const auto [idx, hit] = entry.idxs;
        const bool is_decay = hit == cDecayTrial;
        const bool is_exit = hit == cExitTrial;
        if(!particles[idx].IsPropagating()) continue;
        if(!is_decay && !is_exit && particles[hit].Status() != ParticleStatus::background)
            continue;

        SyncParticles(particles, entry.time);
        currentTime = entry.time;

        if(is_exit) {
            FinishTrack(particles, idx);
            continue;
        }

        const size_t nbefore = particles.size();
        if(is_decay) {
            if(!DecayNow(event, idx)) {
                ScheduleDecay(particles[idx], idx);
                continue;
            }
        } else {
            auto &track = m_tracks.at(idx);
            double energy_before = 0;
            if(track.potential) {
                // Collision kinematics with the momentum on the energy shell at this point
                const auto p = track.path.ShellMomentum(currentTime);
                const double mass = track.path.Mass();
                particles[idx].SetMomentum(
                    {std::sqrt(p.P2() + mass * mass), p.Px(), p.Py(), p.Pz()});
                const auto pot = CurvedPotential();
                const auto &spec = particles[hit];
                energy_before = track.path.Energy() +
                                flow::Hamiltonian(spec.Position(), spec.Momentum().Vec3(), *pot,
                                                  spec.Info().Mass());
            }
            FinalizeMomentum(event, particles, idx, hit);
            if(particles[idx].Status() != ParticleStatus::interacted) {
                ++m_diag.pauli_blocked;
                const auto key = std::make_pair(idx, hit);
                if(m_algorithm == Algorithm::Hazard) {
                    if(PauliMode() == CurvedOptions::Pauli::Thinning) {
                        if(m_kernel.HardSphere()) {
                            // Infinite hazard inside the sphere: retry at once (bounded)
                            if(++m_atom_attempts[key] < 1000)
                                m_time_steps.push({currentTime, {idx, hit}});
                        } else {
                            SchedulePair(event, idx, hit, currentTime);
                        }
                    }
                } else {
                    auto it = m_atom_attempts.find(key);
                    if(it != m_atom_attempts.end() && it->second > 0) {
                        --it->second;
                        m_time_steps.push({currentTime, {idx, hit}});
                    } else {
                        SchedulePair(event, idx, hit, currentTime + 1e-9);
                    }
                }
                continue;
            }
            ++m_diag.collisions;
            if(track.potential) {
                const auto pot = CurvedPotential();
                double energy_after = 0;
                for(size_t i = nbefore; i < particles.size(); ++i) {
                    const auto &out = particles[i];
                    energy_after +=
                        out.Info().IsNucleon()
                            ? flow::Hamiltonian(out.Position(), out.Momentum().Vec3(), *pot,
                                                out.Info().Mass())
                            : out.Momentum().E();
                }
                m_diag.energy_mismatch += std::abs(energy_after - energy_before);
            }
        }
        for(size_t i = nbefore; i < particles.size(); ++i) {
            if(StartTrack(event, i, with_potential)) ScheduleCurvedPairs(event, i);
        }
    }

    // Particles that never reach the exit radius (no finite exit time)
    for(auto &particle : particles) {
        if(particle.IsPropagating()) particle.Status() = ParticleStatus::captured;
    }

    Validate(event);
    Reset();
    m_tracks.clear();
    m_pair_sigma.clear();
}

std::vector<double> Cascade::CumulativeHazard(Track &track, const ThreeVector &x, double sigma,
                                             double start, const std::vector<double> &times) {
    // H(t) = int_start^t Gamma ds at the requested times, by quadrature on the dense path,
    // independent of the sampler's majorants: 8-point Gauss-Legendre on pieces far from the
    // spectator, adaptive Simpson (relative tolerance 1e-10) on pieces that pass close to it
    auto &path = track.path;
    const double texit = path.ExitTime();
    auto f = [&](double t) { return PairHazard(track, x, sigma, t); };
    static constexpr std::array<double, 8> gx{-0.9602898564975363, -0.7966664774136267,
                                              -0.5255324099163290, -0.1834346424956498,
                                              0.1834346424956498,  0.5255324099163290,
                                              0.7966664774136267,  0.9602898564975363};
    static constexpr std::array<double, 8> gw{0.1012285362903763, 0.2223810344533745,
                                              0.3137066278286366, 0.3626837833783620,
                                              0.3626837833783620, 0.3137066278286366,
                                              0.2223810344533745, 0.1012285362903763};
    std::function<double(double, double, double, double, double, double, int)> simpson;
    simpson = [&](double a, double b, double fa, double fm, double fb, double whole,
                  int depth) -> double {
        const double m = 0.5 * (a + b);
        const double flm = f(0.5 * (a + m)), frm = f(0.5 * (m + b));
        const double left = (m - a) / 6 * (fa + 4 * flm + fm);
        const double right = (b - m) / 6 * (fm + 4 * frm + fb);
        if(depth > 50 || std::abs(left + right - whole) <= 1e-10 * std::abs(left + right))
            return left + right + (left + right - whole) / 15;
        return simpson(a, m, fa, flm, fm, left, depth + 1) +
               simpson(m, b, fm, frm, fb, right, depth + 1);
    };
    auto piece = [&](double a, double b) {
        if(b <= a) return 0.0;
        const auto bound = path.Bounds(a, b);
        const double dlow = SegmentDistance(x, bound.a, bound.b) - bound.delta;
        const double length = bound.vmax * (b - a);
        if(dlow > 3 * length) {
            double sum = 0;
            for(size_t i = 0; i < gx.size(); ++i)
                sum += gw[i] * f(0.5 * (a + b) + 0.5 * (b - a) * gx[i]);
            return 0.5 * (b - a) * sum;
        }
        const double fa = f(a), fb = f(b), fm = f(0.5 * (a + b));
        if(!std::isfinite(fa) || !std::isfinite(fb) || !std::isfinite(fm)) return cInf;
        return simpson(a, b, fa, fm, fb, (b - a) / 6 * (fa + 4 * fm + fb), 0);
    };
    std::vector<double> result(times.size(), 0);
    std::vector<size_t> order(times.size());
    for(size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](size_t i, size_t j) { return times[i] < times[j]; });
    double total = 0, t = start;
    for(const size_t i : order) {
        const double target = std::min(times[i], texit);
        while(t < target) {
            const size_t k = path.Segment(t);
            const double tend = std::min(path.SegmentEnd(k), target);
            total += piece(t, tend);
            t = tend;
        }
        result[i] = target > start ? total : 0;
    }
    return result;
}

std::vector<double> Cascade::SurvivalProfile(Event &event, Nucleus *nucleus, size_t idx,
                                             const std::vector<double> &times, size_t npauli) {
    m_nucleus = nucleus;
    currentTime = 0;
    // No state from a previous event: cross sections are cached per (idx, hit) pair
    m_tracks.clear();
    m_pair_sigma.clear();
    std::vector<double> survival(times.size(), 1);
    const bool with_potential = m_potential_prop && m_nucleus->GetPotential() != nullptr &&
                                (m_algorithm == Algorithm::CurvedVeto ||
                                 m_algorithm == Algorithm::Hazard);
    auto made = MakeTrack(event, idx, with_potential, 0);
    if(!made) return std::vector<double>(times.size(), 0);
    m_tracks[idx] = std::move(*made);
    auto &track = m_tracks.at(idx);
    auto &path = track.path;
    const double texit = path.ExitTime();
    Particles &particles = event.Hadrons();
    const bool pauli = npauli > 0 && m_nucleus->ID() != PID::hydrogen() &&
                       m_nucleus->ID() != PID::free_neutron();
    const bool thinning = PauliMode() == CurvedOptions::Pauli::Thinning;

    // Pauli-allowed fraction of pair j with the propagating particle moved to time t
    auto allowed = [&](size_t j, double t) {
        if(!pauli) return 1.0;
        Particle moved = particles[idx];
        moved.SetPosition(path.Position(t));
        const auto p = track.potential ? path.ShellMomentum(t) : moved.Momentum().Vec3();
        const double mass = path.Mass();
        moved.SetMomentum({std::sqrt(p.P2() + mass * mass), p.Px(), p.Py(), p.Pz()});
        auto modes = m_interactions.CrossSection(event, idx, j);
        size_t nallowed = 0;
        for(size_t s = 0; s < npauli; ++s) {
            auto mode = m_interactions.SelectChannel(modes, Random::Instance().Uniform(0.0, 1.0));
            auto outs =
                m_interactions.GenerateMomentum(moved, particles[j], mode.particles, Random::Instance());
            bool blocked = false;
            for(const auto &out : outs) blocked |= PauliBlocking(out);
            nallowed += !blocked;
        }
        return static_cast<double>(nallowed) / static_cast<double>(npauli);
    };

    for(size_t j = 0; j < particles.size(); ++j) {
        if(particles[j].Status() != ParticleStatus::background) continue;
        const double sigma = PairSigma(event, idx, j);
        if(!(sigma > 0)) continue;
        const ThreeVector x = particles[j].Position();
        m_event = &event;
        m_chan_idx = idx;
        m_chan_hit = j;
        m_local_sigma = m_curved.sigma == CurvedOptions::Sigma::Local &&
                        (m_algorithm == Algorithm::Hazard || m_algorithm == Algorithm::CurvedVeto);
        switch(m_algorithm) {
        case Algorithm::Base:
        case Algorithm::MFP:
        case Algorithm::Veto:
        case Algorithm::CurvedVeto: {
            for(const auto &[ta, tb, dlow] : CurvedMinima(track, x, 0, texit)) {
                (void)dlow;
                double tmin = ta;
                if(tb > ta) {
                    double lo = ta, hi = tb;
                    for(size_t it = 0; it < 100; ++it) {
                        const double mid = 0.5 * (lo + hi);
                        ((path.Position(mid) - x).Dot(path.Velocity(mid)) < 0 ? lo : hi) = mid;
                    }
                    tmin = 0.5 * (lo + hi);
                }
                if(tmin < track.t_formed) continue;
                const double sig = m_local_sigma ? LocalSigma(tmin) : sigma;
                const double prob = std::min(1.0, m_kernel.Probability((path.Position(tmin) - x).P2(), sig));
                const double a = prob > 0 ? allowed(j, tmin) : 1;
                const double keep = thinning && m_algorithm == Algorithm::CurvedVeto
                                        ? std::pow(1 - prob, a)
                                        : 1 - prob * a;
                for(size_t k = 0; k < times.size(); ++k)
                    if(tmin <= times[k]) survival[k] *= keep;
            }
            break;
        }
        case Algorithm::Continuous: {
            // Box hazard of SmoothPairTime on the straight line (single attempt)
            const auto v = path.Velocity(0);
            const double speed = v.P();
            const double tc = (x - path.Position(0)).Dot(v) / (speed * speed);
            const auto sep = x - path.Position(0);
            const double along = sep.Dot(v) / speed;
            const double b2 = sep.P2() - along * along;
            const double prob = m_kernel.Probability(b2, sigma);
            if(!(prob > 0)) break;
            const double hw = std::sqrt(3 * sigma / (2 * M_PI)) / speed;
            const double lower = std::max(tc - hw, track.t_formed), upper = tc + hw;
            if(lower >= upper) break;
            const double rate = prob >= 1 ? cInf : -std::log1p(-prob) / (2 * hw);
            const double a = allowed(j, std::clamp(tc, lower, upper));
            for(size_t k = 0; k < times.size(); ++k) {
                const double t = std::min(times[k], texit);
                if(t <= lower) continue;
                const double fired = rate == cInf ? 1 : -std::expm1(-rate * (std::min(t, upper) - lower));
                survival[k] *= 1 - a * fired;
            }
            break;
        }
        case Algorithm::Hazard: {
            if(m_kernel.HardSphere()) {
                auto contact = NextContactTime(idx, x, sigma, track.t_formed);
                if(!contact) break;
                const double a = allowed(j, *contact);
                for(size_t k = 0; k < times.size(); ++k)
                    if(*contact <= times[k]) survival[k] *= thinning ? (a > 0 ? 0 : 1) : 1 - a;
                break;
            }
            // Pauli fraction at the closest point of the path (constant along the pair)
            double tclose = track.t_formed, best = cInf;
            for(double t = track.t_formed; t < texit; t += path.Step()) {
                const double d = (path.Position(t) - x).P();
                if(d < best) {
                    best = d;
                    tclose = t;
                }
            }
            const double a = allowed(j, tclose);
            const auto hz = CumulativeHazard(track, x, sigma, track.t_formed, times);
            for(size_t k = 0; k < times.size(); ++k)
                survival[k] *= thinning ? std::exp(-a * hz[k]) : 1 + a * std::expm1(-hz[k]);
            break;
        }
        }
    }
    m_local_sigma = false;
    m_tracks.clear();
    m_pair_sigma.clear();
    return survival;
}
