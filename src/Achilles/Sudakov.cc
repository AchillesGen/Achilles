// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/Sudakov.hh"
#include "Achilles/Constants.hh"
#include "Achilles/Event.hh"
#include "Achilles/EventHistory.hh"
#include "Achilles/InteractionHandler.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/Random.hh"

#include "spdlog/spdlog.h"
#include "yaml-cpp/yaml.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace achilles::sudakov {

// ==========================================================================================
// Utilities
// ==========================================================================================
double Ein(double x) {
    // Ein(x) = sum_{k>=1} (-1)^{k+1} x^k / (k k!) for small x; E1(x) + ln x + gamma otherwise
    if(x <= 0) return 0;
    if(x < 2.0) {
        double term = x, sum = x;
        for(int k = 2; k < 60; ++k) {
            term *= -x / k;
            const double add = term / k;
            sum += add;
            if(std::abs(add) < 1e-17 * std::abs(sum)) break;
        }
        return sum;
    }
    return std::expint(-x) * -1.0 + std::log(x) + 0.57721566490153286061;
}

namespace {
double SafeU(double u) {
    return u > 0 ? u : 1e-300;
}

double MollerVelocity(const FourVector &p1, const FourVector &p2) {
    const double dot = p1 * p2;
    const double m1 = p1.M2(), m2 = p2.M2();
    const double arg = dot * dot - m1 * m2;
    return arg > 0 ? std::sqrt(arg) / (p1.E() * p2.E()) : 0.0;
}

// Times at which |r0 + beta (t - tref)| = R. Returns (t_in, t_out) or nothing.
bool SphereCrossing(const ThreeVector &r0, const ThreeVector &beta, double tref, double R,
                    double &tin, double &tout) {
    const double a = beta.Magnitude2();
    const double b = r0.Dot(beta);
    const double c = r0.Magnitude2() - R * R;
    const double disc = b * b - a * c;
    if(a <= 0 || disc <= 0) return false;
    const double sq = std::sqrt(disc);
    tin = tref + (-b - sq) / a;
    tout = tref + (-b + sq) / a;
    return true;
}
} // namespace

// ==========================================================================================
// Majorant
// ==========================================================================================
void Majorant::Add(double t0, double t1, double rate) {
    if(!(t1 > t0) || !(rate > 0)) return;
    if(!m_t1.empty() && t0 < m_t1.back() - 1e-12)
        throw std::logic_error("Majorant::Add: segments must be added in time order");
    m_t0.push_back(t0);
    m_t1.push_back(t1);
    m_rate.push_back(rate);
}

void Majorant::SortAtoms() {
    std::sort(m_atoms.begin(), m_atoms.end(),
              [](const Atom &a, const Atom &b) { return a.time < b.time; });
}

double Majorant::NextContinuous(double t, double u) const {
    double target = -std::log(SafeU(u));
    for(size_t i = 0; i < m_rate.size(); ++i) {
        if(m_t1[i] <= t) continue;
        const double start = std::max(t, m_t0[i]);
        const double avail = m_rate[i] * (m_t1[i] - start);
        if(avail >= target) return start + target / m_rate[i];
        target -= avail;
    }
    return kInf;
}

const Majorant::Atom *Majorant::NextAtom(double t) const {
    auto it = std::upper_bound(m_atoms.begin(), m_atoms.end(), t,
                               [](double tt, const Atom &a) { return tt < a.time; });
    return it == m_atoms.end() ? nullptr : &*it;
}

double Majorant::RateAt(double t) const {
    auto it = std::upper_bound(m_t0.begin(), m_t0.end(), t);
    if(it == m_t0.begin()) return 0;
    const size_t i = static_cast<size_t>(std::distance(m_t0.begin(), it)) - 1;
    return t < m_t1[i] ? m_rate[i] : 0.0;
}

double Majorant::Integral(double t0, double t1) const {
    double sum = 0;
    for(size_t i = 0; i < m_rate.size(); ++i) {
        const double a = std::max(t0, m_t0[i]), b = std::min(t1, m_t1[i]);
        if(b > a) sum += m_rate[i] * (b - a);
    }
    return sum;
}

double Majorant::End() const {
    double end = m_t1.empty() ? -kInf : m_t1.back();
    if(!m_atoms.empty()) end = std::max(end, m_atoms.back().time);
    return end;
}

// ==========================================================================================
// Occupancy
// ==========================================================================================
Occupancy::Occupancy(Shape shape, double delta) : m_shape(shape), m_delta(delta) {
    Tabulate();
}

void Occupancy::Tabulate() {
    // Phi normalised so that 3 int Phi(x) x^2 dx = 1
    if(m_shape == Shape::Step) {
        m_norm = 1;
        m_xmax = 1;
    } else {
        m_xmax = 1 + 25 * m_delta;
        const size_t n = 4001;
        double sum = 0;
        for(size_t i = 0; i < n; ++i) {
            const double x = m_xmax * static_cast<double>(i) / static_cast<double>(n - 1);
            const double w = (i == 0 || i == n - 1) ? 0.5 : 1.0;
            sum += w * 3 * x * x / (1 + std::exp((x - 1) / m_delta));
        }
        sum *= m_xmax / static_cast<double>(n - 1);
        m_norm = 1 / sum;
        if(m_norm > 1 + 1e-9)
            throw std::runtime_error("Occupancy: FermiDirac normalisation requires Z > 1");
    }
    const size_t n = 2001;
    m_xgrid.resize(n);
    m_cdf.resize(n);
    double acc = 0;
    for(size_t i = 0; i < n; ++i) {
        m_xgrid[i] = m_xmax * static_cast<double>(i) / static_cast<double>(n - 1);
        if(i > 0) {
            const double xm = 0.5 * (m_xgrid[i] + m_xgrid[i - 1]);
            acc += 3 * xm * xm * Phi(xm) * (m_xgrid[i] - m_xgrid[i - 1]);
        }
        m_cdf[i] = acc;
    }
    for(auto &c : m_cdf) c /= acc;
}

void Occupancy::SetDensities(std::function<double(double)> rho_p,
                             std::function<double(double)> rho_n) {
    m_rho_p = std::move(rho_p);
    m_rho_n = std::move(rho_n);
}

double Occupancy::KFermi(double r, PID pid) const {
    const double rho = pid == PID::proton() ? m_rho_p(r) : m_rho_n(r);
    return rho > 0 ? std::cbrt(3 * M_PI * M_PI * rho) * Constant::HBARC : 0.0;
}

double Occupancy::Phi(double x) const {
    if(m_shape == Shape::Step) return x < 1 ? 1.0 : 0.0;
    return m_norm / (1 + std::exp((x - 1) / m_delta));
}

double Occupancy::f(double r, PID pid, double k) const {
    if(pid != PID::proton() && pid != PID::neutron()) return 0;
    const double kf = KFermi(r, pid);
    if(kf <= 0) return 0;
    return Phi(k / kf);
}

double Occupancy::SampleK(double r, PID pid, double u) const {
    auto it = std::lower_bound(m_cdf.begin(), m_cdf.end(), u);
    size_t i = static_cast<size_t>(std::distance(m_cdf.begin(), it));
    if(i == 0) i = 1;
    if(i >= m_cdf.size()) i = m_cdf.size() - 1;
    const double c0 = m_cdf[i - 1], c1 = m_cdf[i];
    const double x =
        m_xgrid[i - 1] + (c1 > c0 ? (u - c0) / (c1 - c0) : 0.0) * (m_xgrid[i] - m_xgrid[i - 1]);
    return x * KFermi(r, pid);
}

// ==========================================================================================
// Engine
// ==========================================================================================
Engine::Engine(std::shared_ptr<RateModel> rate, std::shared_ptr<CollisionModel> coll,
               EngineSettings settings, Uniform01 uniform)
    : m_rate(std::move(rate)), m_coll(std::move(coll)), m_settings(std::move(settings)) {
    m_ctx.uniform = std::move(uniform);
    m_ctx.diag = &m_diag;
}

void Engine::Push(const Entry &e) {
    m_heap.push_back(e);
    std::push_heap(m_heap.begin(), m_heap.end(), std::greater<Entry>{});
}

Engine::Entry Engine::Pop() {
    std::pop_heap(m_heap.begin(), m_heap.end(), std::greater<Entry>{});
    Entry e = m_heap.back();
    m_heap.pop_back();
    return e;
}

void Engine::PropagateTo(Event &event, double t) {
    const double dt = t - m_ctx.time;
    if(dt > 0) {
        for(auto &p : event.Hadrons())
            if(p.IsPropagating()) p.Propagate(dt);
    }
    m_ctx.time = t;
}

void Engine::Schedule(Event &event, size_t idx, double t0, bool rebuild, bool new_line) {
    auto &part = event.Hadrons()[idx];
    // Formation zone: the rate is switched on after the remaining formation time
    double start = t0;
    if(part.InFormationZone() && t0 <= m_ctx.time + 1e-12) start = t0 + part.FormationZone();
    if(rebuild) {
        m_ctx.new_line = new_line;
        m_major[idx] = m_rate->BuildMajorant(event, idx, start, m_ctx);
    }
    const Majorant &M = m_major[idx];
    const double tc = M.NextContinuous(start, m_ctx.uniform());
    const Majorant::Atom *atom = M.NextAtom(start);
    if(atom && atom->time < tc)
        Push({atom->time, idx, m_gen[idx], true, atom->label});
    else if(tc < kInf)
        Push({tc, idx, m_gen[idx], false, 0});
}

Engine::~Engine() {
    // run summary; majorant violations mean a bound (e.g. MediumMax) was too small
    try {
        if(m_diag.events > 0)
            spdlog::info("Sudakov: {} events, {} trials, {} accepted, {} majorant violations "
                         "(max factor {:.3f})",
                         m_diag.events, m_diag.trials, m_diag.accepted, m_diag.violations,
                         m_diag.max_ratio);
    } catch(...) {
    }
}

void Engine::Evolve(Event &event, Nucleus *nucleus) {
    m_ctx.time = 0;
    m_ctx.n_struck = 0;
    m_ctx.collisions.clear();
    m_ctx.hole_centres.clear();
    m_ctx.weights.assign(m_settings.variations.size(), 1.0);
    m_diag.events++;
    m_heap.clear();

    m_coll->BeginEvent(event, nucleus);
    m_rate->BeginEvent(event, nucleus, m_ctx);

    auto &parts = event.Hadrons();
    m_gen.assign(parts.size(), 0);
    m_major.assign(parts.size(), Majorant{});
    for(size_t i = 0; i < parts.size(); ++i) {
        if(parts[i].IsPropagating() && m_coll->Interacts(parts[i]))
            Schedule(event, i, 0, true, true);
    }

    const size_t nvar = m_settings.variations.size();
    size_t ntrials = 0;
    while(!m_heap.empty()) {
        const Entry e = Pop();
        if(e.idx >= m_gen.size() || e.gen != m_gen[e.idx]) continue;
        if(!event.Hadrons()[e.idx].IsPropagating()) continue;
        if(++ntrials > m_settings.max_trials_per_event) {
            m_diag.capped++;
            spdlog::warn("Sudakov: trial cap reached in event {}", m_diag.events);
            break;
        }
        PropagateTo(event, e.time);

        Trial trial;
        trial.proj = e.idx;
        trial.time = e.time;
        trial.atom = e.atom;
        trial.atom_label = e.atom_label;
        trial.point = event.Hadrons()[e.idx].Position();
        trial.majorant_rate = e.atom ? 1.0 : m_major[e.idx].RateAt(e.time);
        bool possible = m_rate->Propose(event, trial, m_ctx);
        trial.n_rate_factors = trial.factors.size();
        if(possible) m_coll->Propose(event, trial, m_ctx);

        // two stages: rate-model factors (a1) and collision-model factors (a2)
        double a1 = possible ? 1.0 : 0.0, a2 = 1.0;
        std::vector<double> av1(nvar, a1), av2(nvar, 1.0);
        for(size_t k = 0; k < trial.factors.size(); ++k) {
            auto &f = trial.factors[k];
            const bool stage1 = k < trial.n_rate_factors;
            if(f.value > 1 + 1e-12 || f.value < -1e-12) {
                m_diag.violations++;
                m_diag.max_ratio = std::max(m_diag.max_ratio, f.value);
                spdlog::debug("Sudakov: factor {} = {} outside [0,1] at t = {}", f.name, f.value,
                              e.time);
                if(m_settings.clip_violations) f.value = std::clamp(f.value, 0.0, 1.0);
            }
            (stage1 ? a1 : a2) *= f.value;
            for(size_t v = 0; v < nvar; ++v) {
                double fv = f.vars.empty() ? f.value : f.vars[v];
                if(fv > 1 + 1e-12) {
                    m_diag.violations++;
                    m_diag.max_ratio = std::max(m_diag.max_ratio, fv);
                }
                (stage1 ? av1[v] : av2[v]) *= std::clamp(fv, 0.0, 1.0);
            }
        }
        m_diag.trials++;
        const bool pass1 = possible && m_ctx.uniform() < a1;
        const bool pass2 = pass1 && m_ctx.uniform() < a2;
        const bool accept = pass2;
        // stage-resolved likelihood ratio (equals the one-stage formula when nothing
        // depends on which stage rejected)
        for(size_t v = 0; v < nvar; ++v) {
            if(accept)
                m_ctx.weights[v] *= (av1[v] * av2[v]) / (a1 * a2);
            else if(!pass1) {
                if(a1 < 1) m_ctx.weights[v] *= (1 - av1[v]) / (1 - a1);
            } else if(a2 < 1)
                m_ctx.weights[v] *= av1[v] * (1 - av2[v]) / (a1 * (1 - a2));
        }

        if(!accept) {
            if(pass1 && m_settings.single_attempt) m_rate->OnCollisionVeto(event, trial, m_ctx);
            m_rate->OnReject(event, trial, m_ctx);
            Schedule(event, e.idx, e.time, false, false);
            continue;
        }

        m_diag.accepted++;
        m_rate->OnAccept(event, trial, m_ctx);
        const auto proj_id = event.Hadrons()[trial.proj].ID();
        const auto targ_id = event.Hadrons()[trial.target].ID();
        auto newparts = m_coll->Commit(event, trial, m_ctx);
        m_ctx.collisions.push_back(
            {e.time, trial.proj, trial.target, trial.point, proj_id, targ_id});

        const size_t n = event.Hadrons().size();
        m_gen.resize(n, 0);
        m_major.resize(n);
        const bool invalid = m_rate->CommitInvalidatesOthers();
        if(event.Hadrons()[e.idx].IsPropagating()) {
            m_gen[e.idx]++;
            Schedule(event, e.idx, e.time, trial.deflects || invalid, trial.deflects);
        }
        for(size_t k : newparts) {
            if(event.Hadrons()[k].IsPropagating() && m_coll->Interacts(event.Hadrons()[k]))
                Schedule(event, k, e.time, true, true);
        }
        if(invalid) {
            for(size_t k = 0; k < n; ++k) {
                if(k == e.idx || !event.Hadrons()[k].IsPropagating() ||
                   !m_coll->Interacts(event.Hadrons()[k]))
                    continue;
                m_gen[k]++;
                Schedule(event, k, e.time, true, false);
            }
        }
    }
}

// ==========================================================================================
// Track A: density backbone
// ==========================================================================================
DensityRate::DensityRate(DensityRateSettings s, std::shared_ptr<CollisionModel> coll,
                         std::shared_ptr<Occupancy> occ)
    : m_s(std::move(s)), m_coll(std::move(coll)), m_occ(std::move(occ)) {}

void DensityRate::SetDensities(std::function<double(double)> rho_p,
                               std::function<double(double)> rho_n, double A, double rmax) {
    m_rho_p = std::move(rho_p);
    m_rho_n = std::move(rho_n);
    m_A = A;
    m_rmax = rmax;
    m_ready = false;
}

void DensityRate::Tabulate() {
    const double rmax0 = m_rmax;
    const double w = m_s.fold_width;
    const double rmax = rmax0 + (w > 0 ? 5 * w : 0);
    const size_t n = 2001;
    m_r.resize(n);
    m_rho_tab.resize(n);
    m_rhop_tab.resize(n);
    m_rhon_tab.resize(n);
    // raw densities truncated at rmax0
    auto rp = [&](double r) { return r < rmax0 ? m_rho_p(r) : 0.0; };
    auto rn = [&](double r) { return r < rmax0 ? m_rho_n(r) : 0.0; };
    for(size_t i = 0; i < n; ++i)
        m_r[i] = rmax * static_cast<double>(i) / static_cast<double>(n - 1);
    if(w <= 0) {
        for(size_t i = 0; i < n; ++i) {
            m_rhop_tab[i] = rp(m_r[i]);
            m_rhon_tab[i] = rn(m_r[i]);
        }
    } else {
        // spherical 3D Gaussian folding:
        // rho_s(r) = 1/(sqrt(2 pi) w r) int r' rho(r') [e^{-(r-r')^2/2w^2} - e^{-(r+r')^2/2w^2}]
        // dr'
        const size_t nq = 3001;
        const double h = rmax0 / static_cast<double>(nq - 1);
        for(size_t i = 0; i < n; ++i) {
            const double r = std::max(m_r[i], 1e-4);
            double sp = 0, sn = 0;
            for(size_t q = 0; q < nq; ++q) {
                const double rr = h * static_cast<double>(q);
                const double wq = (q == 0 || q == nq - 1) ? 0.5 : 1.0;
                const double kern = rr * (std::exp(-(r - rr) * (r - rr) / (2 * w * w)) -
                                          std::exp(-(r + rr) * (r + rr) / (2 * w * w)));
                sp += wq * kern * rp(rr);
                sn += wq * kern * rn(rr);
            }
            const double pre = h / (std::sqrt(2 * M_PI) * w * r);
            m_rhop_tab[i] = pre * sp;
            m_rhon_tab[i] = pre * sn;
        }
    }
    for(size_t i = 0; i < n; ++i) m_rho_tab[i] = m_rhop_tab[i] + m_rhon_tab[i];

    const size_t ns = std::max<size_t>(1, m_s.shells);
    m_shell_edges.resize(ns + 1);
    m_shell_max.assign(ns, 0);
    for(size_t k = 0; k <= ns; ++k)
        m_shell_edges[k] = rmax * static_cast<double>(k) / static_cast<double>(ns);
    for(size_t i = 0; i < n; ++i) {
        size_t k = std::min(ns - 1, static_cast<size_t>(m_r[i] / rmax * static_cast<double>(ns)));
        m_shell_max[k] = std::max(m_shell_max[k], m_rho_tab[i]);
        // grid points on an edge also bound the neighbouring shell
        if(k > 0 && std::abs(m_r[i] - m_shell_edges[k]) < 1e-12)
            m_shell_max[k - 1] = std::max(m_shell_max[k - 1], m_rho_tab[i]);
    }
    for(auto &m : m_shell_max) m *= 1.01; // safety for interpolation between grid points
    if(m_s.rho_hat_override > 0)
        for(auto &m : m_shell_max) m = m_s.rho_hat_override;
    m_ready = true;
}

double DensityRate::RhoTotal(double r) const {
    if(r >= m_r.back()) return 0;
    const double x = r / m_r.back() * static_cast<double>(m_r.size() - 1);
    const size_t i = static_cast<size_t>(x);
    const double f = x - static_cast<double>(i);
    return (1 - f) * m_rho_tab[i] + f * m_rho_tab[i + 1];
}

double DensityRate::HoleG(double d) const {
    if(m_s.hole_d.empty()) return 1;
    if(d >= m_s.hole_d.back()) return 1;
    auto it = std::upper_bound(m_s.hole_d.begin(), m_s.hole_d.end(), d);
    if(it == m_s.hole_d.begin()) return m_s.hole_g.front();
    const size_t i = static_cast<size_t>(std::distance(m_s.hole_d.begin(), it));
    const double f = (d - m_s.hole_d[i - 1]) / (m_s.hole_d[i] - m_s.hole_d[i - 1]);
    return (1 - f) * m_s.hole_g[i - 1] + f * m_s.hole_g[i];
}

void DensityRate::BeginEvent(Event &event, Nucleus *nuc, Context &ctx) {
    m_nuc = nuc;
    if(!m_rho_p && nuc) {
        m_rho_p = [nuc](double r) { return nuc->ProtonRho(r); };
        m_rho_n = [nuc](double r) { return nuc->NeutronRho(r); };
        m_A = static_cast<double>(nuc->NNucleons());
        if(m_rmax <= 0) m_rmax = m_s.rmax > 0 ? m_s.rmax : nuc->Radius();
        m_ready = false;
    }
    if(m_rmax <= 0 && m_s.rmax > 0) m_rmax = m_s.rmax;
    if(!m_ready) {
        if(!m_rho_p) throw std::runtime_error("SudakovA: no density supplied");
        Tabulate();
        if(!m_occ->HasDensities()) m_occ->SetDensities(m_rho_p, m_rho_n);
        // largest Fermi momentum of the occupancy used to sample targets (bounds v_rel, sigma)
        m_kfmax = 0;
        for(size_t i = 0; i < m_r.size(); ++i)
            m_kfmax = std::max({m_kfmax, m_occ->KFermi(m_r[i], PID::proton()),
                                m_occ->KFermi(m_r[i], PID::neutron())});
    }
    m_lines.assign(event.Hadrons().size(), Line{});
    if(m_s.vertex_hole && !m_s.hole_d.empty()) {
        for(const auto &p : event.Hadrons())
            if(p.Status() == ParticleStatus::internal_test)
                ctx.hole_centres.push_back(p.Position());
    }
}

double DensityRate::Q(const Line &L, double t) const {
    if(!m_s.finite_A) return 0;
    // (sigma_ref / A) int_{t0}^{t} rho_s(|r(t')|) v dt'
    const double tt = std::max(t, L.t0);
    const size_t n = 64;
    const double h = (tt - L.t0) / static_cast<double>(n);
    if(h <= 0) return 0;
    double sum = 0;
    for(size_t i = 0; i <= n; ++i) {
        const double ti = L.t0 + h * static_cast<double>(i);
        const double w = (i == 0 || i == n) ? 1 : (i % 2 ? 4 : 2);
        sum += w * RhoTotal((L.r0 + L.beta * (ti - L.t0)).Magnitude());
    }
    return L.sigma_ref * kMbToFm2 * sum * h / 3 * L.beta.Magnitude() / m_A;
}

Majorant DensityRate::BuildMajorant(Event &event, size_t idx, double t0, Context &ctx) {
    if(m_lines.size() < event.Hadrons().size()) m_lines.resize(event.Hadrons().size());
    const Particle p = event.Hadrons()[idx]; // copy: SigmaBound may append temporaries
    const ThreeVector beta = p.Beta();
    const double v = beta.Magnitude();
    Majorant M;
    if(v <= 0) return M;
    Line &L = m_lines[idx];
    // a new straight line starts at t0 (the particle sits at ctx.time); otherwise keep the
    // line origin so that the finite-A void integral Q continues
    if(ctx.new_line || L.t0 > t0) {
        L.r0 = p.Position() + beta * (t0 - ctx.time);
        L.t0 = t0;
    }
    L.beta = beta;
    const double kmax = m_occ->GetShape() == Occupancy::Shape::Step
                            ? 1.0001 * m_kfmax
                            : m_s.kmax_factor * m_kfmax * m_occ->XMax();
    double smax = 1;
    for(double s : m_s.sigma_scale_vars) smax = std::max(smax, s);
    L.sigma_hat = smax * m_coll->SigmaBound(event, p, kmax);
    const double mN = Constant::mN;
    L.vhat = v + kmax / std::sqrt(kmax * kmax + mN * mN);
    L.sigma_ref = m_coll->SigmaAtRest(event, p);

    double tin, tout;
    const double rmax = m_r.back();
    const ThreeVector rt0 = p.Position() + beta * (t0 - ctx.time);
    if(!SphereCrossing(rt0, beta, t0, rmax, tin, tout) || tout <= t0) return M;
    tin = std::max(tin, t0);

    // bounds of history-dependent factors
    double gfac = 1;
    if(!m_s.hole_d.empty())
        gfac = std::pow(m_s.hole_gmax, static_cast<double>(ctx.hole_centres.size()));
    double void_bound = 1;
    if(m_s.finite_A) {
        // the IPM void factor 1/(1-Q) assumes a linear per-nucleon probability Q = sigma T / A;
        // it is capped at Q = kQcap (model limit, counted in diagnostics)
        const double qexit = std::min(Q(L, tout), kQcap);
        L.qexit = qexit;
        void_bound = 1 / (1 - qexit);
    }
    const double pre = m_s.headroom * L.sigma_hat * kMbToFm2 * L.vhat * gfac * void_bound;

    // breakpoints at the shell edges
    std::vector<double> times{tin, tout};
    for(size_t k = 1; k + 1 < m_shell_edges.size(); ++k) {
        double a, b;
        if(SphereCrossing(rt0, beta, t0, m_shell_edges[k], a, b)) {
            if(a > tin && a < tout) times.push_back(a);
            if(b > tin && b < tout) times.push_back(b);
        }
    }
    std::sort(times.begin(), times.end());
    const size_t ns = m_shell_max.size();
    for(size_t i = 0; i + 1 < times.size(); ++i) {
        const double ta = times[i], tb = times[i + 1];
        if(tb - ta <= 1e-12) continue;
        const double rmid = (rt0 + beta * (0.5 * (ta + tb) - t0)).Magnitude();
        const size_t k =
            std::min(ns - 1, static_cast<size_t>(rmid / rmax * static_cast<double>(ns)));
        M.Add(ta, tb, pre * m_shell_max[k]);
    }
    return M;
}

bool DensityRate::Propose(Event &event, Trial &trial, Context &ctx) {
    const auto &L = m_lines[trial.proj];
    const double r = trial.point.Magnitude();
    const double rho = RhoTotal(r);
    if(rho <= 0) return false;
    // species mark ~ rho_tau / rho (folded tables)
    const double x = r / m_r.back() * static_cast<double>(m_r.size() - 1);
    const size_t i = std::min(static_cast<size_t>(x), m_r.size() - 2);
    const double fr = x - static_cast<double>(i);
    const double rp = (1 - fr) * m_rhop_tab[i] + fr * m_rhop_tab[i + 1];
    const PID pid = ctx.uniform() < rp / rho ? PID::proton() : PID::neutron();
    // momentum mark ~ f(r, k) (unfolded local density at the point)
    const double k = m_occ->SampleK(r, pid, ctx.uniform());
    const double cth = 2 * ctx.uniform() - 1, phi = 2 * M_PI * ctx.uniform();
    const double sth = std::sqrt(std::max(0.0, 1 - cth * cth));
    const double m = ParticleInfo(pid).Mass();
    FourVector mom{std::sqrt(m * m + k * k), k * sth * std::cos(phi), k * sth * std::sin(phi),
                   k * cth};
    event.Hadrons().emplace_back(pid, mom, trial.point, ParticleStatus::background);
    trial.target = event.Hadrons().size() - 1;
    trial.created_target = true;

    const auto &proj = event.Hadrons()[trial.proj];
    const double sigma = m_coll->SigmaTot(event, trial.proj, trial.target) *
                         m_coll->Medium(proj, event.Hadrons()[trial.target]);
    const double vrel = MollerVelocity(proj.Momentum(), mom);
    double rate = sigma * kMbToFm2 * vrel * rho;
    if(!m_s.hole_d.empty()) {
        for(const auto &c : ctx.hole_centres) rate *= HoleG((trial.point - c).Magnitude());
    }
    if(m_s.finite_A) {
        const double depl = std::max(0.0, (m_A - static_cast<double>(ctx.n_struck)) / m_A);
        rate *= depl / (1 - std::min(Q(L, trial.time), kQcap));
    }
    const double hat = trial.majorant_rate;
    Factor f{"rate", hat > 0 ? rate / hat : 0.0, {}};
    if(!m_s.sigma_scale_vars.empty()) {
        f.vars.resize(m_s.sigma_scale_vars.size());
        for(size_t v = 0; v < f.vars.size(); ++v) f.vars[v] = f.value * m_s.sigma_scale_vars[v];
    }
    trial.factors.push_back(std::move(f));
    return true;
}

void DensityRate::OnAccept(Event &, const Trial &trial, Context &ctx) {
    ctx.n_struck++;
    if(!m_s.hole_d.empty()) ctx.hole_centres.push_back(trial.point);
}

void DensityRate::OnReject(Event &event, const Trial &trial, Context &) {
    if(trial.created_target && trial.target == event.Hadrons().size() - 1)
        event.Hadrons().pop_back();
}

// ==========================================================================================
// Track B: configuration backbone
// ==========================================================================================
ConfigurationRate::ConfigurationRate(ConfigurationRateSettings s,
                                     std::shared_ptr<CollisionModel> coll)
    : m_s(std::move(s)), m_coll(std::move(coll)) {}

double ConfigurationRate::Width2(double sigma) const {
    if(m_s.law == ConfigurationRateSettings::Law::Kappa1) return sigma / (2 * M_PI);
    return sigma / (2 * M_PI * Ein(m_s.x0));
}

double ConfigurationRate::PairX(double b2, double sigma) const {
    if(sigma <= 0) return 0;
    if(m_s.law == ConfigurationRateSettings::Law::Kappa1) {
        const double p = std::exp(-M_PI * b2 / sigma);
        return p >= 1 - 1e-15 ? 34.5 : -std::log1p(-p);
    }
    return m_s.x0 * std::exp(-b2 / (2 * Width2(sigma)));
}

double ConfigurationRate::PairProb(double b2, double sigma) const {
    return -std::expm1(-PairX(b2, sigma));
}

namespace {
// truncated, renormalised Gaussian in time
double GaussPhi(double dt, double st, double nsig) {
    if(std::abs(dt) > nsig * st) return 0;
    const double norm = std::erf(nsig / std::sqrt(2.0));
    return std::exp(-dt * dt / (2 * st * st)) / (std::sqrt(2 * M_PI) * st * norm);
}
} // namespace

double ConfigurationRate::PairHazard(const Pair &P, double sigma, double t) const {
    const double X = PairX(P.b2, sigma);
    const double s = std::sqrt(Width2(sigma));
    switch(m_s.profile) {
    case ConfigurationRateSettings::Profile::Box: {
        const double lt = std::sqrt(3.0) * s / P.v;
        return std::abs(t - P.tc) < lt ? X / (2 * lt) : 0.0;
    }
    case ConfigurationRateSettings::Profile::Gauss:
        return X * GaussPhi(t - P.tc, s / P.v, m_s.gauss_nsigma);
    default:
        return 0;
    }
}

void ConfigurationRate::PairBoundPieces(const Pair &P, double t0,
                                        std::vector<std::array<double, 3>> &out) const {
    // sum over sigma variants (nominal + variations) of an upper envelope of each hazard
    for(double sigma : P.sigmas) {
        const double X = PairX(P.b2, sigma);
        const double s = std::sqrt(Width2(sigma));
        if(m_s.profile == ConfigurationRateSettings::Profile::Box) {
            const double lt = std::sqrt(3.0) * s / P.v;
            const double a = std::max(t0, P.tc - lt), b = P.tc + lt;
            if(b > a) out.push_back({a, b, m_s.headroom * X / (2 * lt)});
        } else {
            const double st = s / P.v, n = m_s.gauss_nsigma;
            const size_t K = m_s.gauss_bins;
            const double lo = P.tc - n * st, w = 2 * n * st / static_cast<double>(K);
            for(size_t k = 0; k < K; ++k) {
                const double a0 = lo + w * static_cast<double>(k), b0 = a0 + w;
                const double a = std::max(t0, a0);
                if(b0 <= a) continue;
                // max of the Gaussian on [a0, b0]: at tc if inside, else at the nearer edge
                const double tmax = std::clamp(P.tc, a0, b0);
                out.push_back({a, b0, m_s.headroom * X * GaussPhi(tmax - P.tc, st, n + 1e-9)});
            }
        }
    }
}

double ConfigurationRate::PairBound(const Pair &P, double t) const {
    std::vector<std::array<double, 3>> pieces;
    PairBoundPieces(P, -kInf, pieces);
    double sum = 0;
    for(const auto &pc : pieces)
        if(t >= pc[0] && t < pc[1]) sum += pc[2];
    return sum;
}

Majorant ConfigurationRate::BuildMajorant(Event &event, size_t idx, double t0, Context &ctx) {
    if(m_cache.size() < event.Hadrons().size()) m_cache.resize(event.Hadrons().size());
    auto &cache = m_cache[idx];
    cache.pairs.clear();
    const auto &p = event.Hadrons()[idx];
    const ThreeVector beta = p.Beta();
    const double v = beta.Magnitude();
    Majorant M;
    if(v <= 0) return M;
    const ThreeVector r0 = p.Position(); // at ctx.time
    std::vector<std::array<double, 3>> pieces;
    for(size_t j = 0; j < event.Hadrons().size(); ++j) {
        const auto &q = event.Hadrons()[j];
        if(j == idx || !q.IsBackground() || !m_coll->Interacts(q)) continue;
        const ThreeVector d = q.Position() - r0;
        const double along = d.Dot(beta) / v;
        const double b2 = std::max(0.0, d.Magnitude2() - along * along);
        Pair P{j, ctx.time + along / v, b2, v, {}, 0, 0};
        double sig = m_coll->SigmaTot(event, idx, j) * kMbToFm2;
        if(sig > 0) {
            Particle pc = p; // projectile at closest approach to q
            pc.SetPosition(r0 + beta * (along / v));
            sig *= m_coll->Medium(pc, q);
        }
        P.sigmas.push_back(sig);
        for(double s : m_s.sigma_scale_vars) P.sigmas.push_back(s * sig);
        double xmax = 0;
        for(double s : P.sigmas) xmax = std::max(xmax, PairX(b2, s));
        if(xmax < 1e-9) continue;
        if(m_s.profile == ConfigurationRateSettings::Profile::Atom) {
            if(P.tc > t0) {
                M.AddAtom(P.tc, cache.pairs.size());
                cache.pairs.push_back(std::move(P));
            }
            continue;
        }
        const size_t before = pieces.size();
        PairBoundPieces(P, t0, pieces);
        if(pieces.size() > before) cache.pairs.push_back(std::move(P));
    }
    M.SortAtoms();
    if(pieces.empty()) return M;
    // sweep: sum overlapping pieces into non-overlapping segments
    std::vector<std::pair<double, double>> ev;
    ev.reserve(2 * pieces.size());
    for(const auto &pc : pieces) {
        ev.push_back({pc[0], pc[2]});
        ev.push_back({pc[1], -pc[2]});
    }
    std::sort(ev.begin(), ev.end(), [](const auto &a, const auto &b) { return a.first < b.first; });
    double rate = 0;
    for(size_t i = 0; i < ev.size(); ++i) {
        rate += ev[i].second;
        if(i + 1 < ev.size() && ev[i + 1].first > ev[i].first && rate > 1e-14)
            M.Add(ev[i].first, ev[i + 1].first, rate);
    }
    return M;
}

bool ConfigurationRate::Propose(Event &event, Trial &trial, Context &ctx) {
    auto &cache = m_cache[trial.proj];
    if(trial.atom) {
        const auto &P = cache.pairs[trial.atom_label];
        if(!event.Hadrons()[P.j].IsBackground()) return false;
        trial.target = P.j;
        trial.pair = trial.atom_label;
        Factor f{"pair", -std::expm1(-PairX(P.b2, P.sigmas[0])), {}};
        for(size_t v = 1; v < P.sigmas.size(); ++v)
            f.vars.push_back(-std::expm1(-PairX(P.b2, P.sigmas[v])));
        trial.factors.push_back(std::move(f));
        return true;
    }
    // choose the pair with probability hhat_j(t) / sum hhat
    std::vector<double> hb(cache.pairs.size());
    double H = 0;
    for(size_t i = 0; i < cache.pairs.size(); ++i) {
        hb[i] = PairBound(cache.pairs[i], trial.time);
        H += hb[i];
    }
    if(H <= 0) return false;
    double u = ctx.uniform() * H;
    size_t sel = 0;
    for(; sel + 1 < hb.size(); ++sel) {
        if(u < hb[sel]) break;
        u -= hb[sel];
    }
    const auto &P = cache.pairs[sel];
    if(!event.Hadrons()[P.j].IsBackground() || P.blocked) return false;
    trial.target = P.j;
    trial.pair = sel;
    Factor f{"pair", PairHazard(P, P.sigmas[0], trial.time) / hb[sel], {}};
    for(size_t v = 1; v < P.sigmas.size(); ++v)
        f.vars.push_back(PairHazard(P, P.sigmas[v], trial.time) / hb[sel]);
    trial.factors.push_back(std::move(f));
    return true;
}

void ConfigurationRate::OnCollisionVeto(Event &, const Trial &trial, Context &) {
    auto &cache = m_cache[trial.proj];
    if(trial.pair < cache.pairs.size()) cache.pairs[trial.pair].blocked = true;
}

// ==========================================================================================
// Collision models
// ==========================================================================================
void ForwardCollision::Propose(Event &, Trial &trial, Context &) {
    trial.deflects = false;
    if(m_a < 1 || !m_avars.empty()) trial.factors.push_back({"pauli", m_a, m_avars});
}

std::vector<size_t> ForwardCollision::Commit(Event &event, Trial &trial, Context &) {
    event.Hadrons()[trial.target].Status() = ParticleStatus::interacted;
    return {};
}

void AchillesCollision::BeginEvent(Event &, Nucleus *nuc) {
    if(nuc && !m_occ->HasDensities()) {
        m_occ->SetDensities([nuc](double r) { return nuc->ProtonRho(r); },
                            [nuc](double r) { return nuc->NeutronRho(r); });
    }
    if(nuc && m_s.coherent && m_s.veto_q.empty() && !HasFormFactors()) {
        SetFormFactors([nuc](double r) { return nuc->ProtonRho(r); },
                       [nuc](double r) { return nuc->NeutronRho(r); });
    }
}

void AchillesCollision::SetFormFactors(const std::function<double(double)> &rho_p,
                                       const std::function<double(double)> &rho_n,
                                       double rmax) {
    // F(q) = int rho(r) j0(q r) r^2 dr / int rho(r) r^2 dr on a q grid of kFFdq MeV
    const size_t nr = 3001, nq = 2501;
    const double dr = rmax / static_cast<double>(nr - 1);
    for(int species = 0; species < 2; ++species) {
        const auto &rho = species == 0 ? rho_p : rho_n;
        std::vector<double> w(nr);
        double norm = 0;
        for(size_t i = 0; i < nr; ++i) {
            const double r = dr * static_cast<double>(i);
            const double simpson = (i == 0 || i == nr - 1) ? 1.0 : (i % 2 ? 4.0 : 2.0);
            w[i] = simpson * rho(r) * r * r;
            norm += w[i];
        }
        auto &ff = species == 0 ? m_ff_p : m_ff_n;
        ff.assign(nq, 0.0);
        if(norm <= 0) continue;
        for(size_t iq = 0; iq < nq; ++iq) {
            const double q = kFFdq * static_cast<double>(iq) / Constant::HBARC; // fm^-1
            double sum = 0;
            for(size_t i = 0; i < nr; ++i) {
                const double x = q * dr * static_cast<double>(i);
                sum += w[i] * (x < 1e-6 ? 1 - x * x / 6 : std::sin(x) / x);
            }
            ff[iq] = sum / norm;
        }
    }
}

double AchillesCollision::CoherentVetoProb(PID pid, double q) const {
    if(m_s.veto_q.empty()) {
        const double F = FormFactor(pid, q);
        return F * F;
    }
    const auto &xq = m_s.veto_q;
    const auto &yv = m_s.veto_v;
    if(q <= xq.front()) return std::clamp(yv.front(), 0.0, 1.0);
    if(q >= xq.back()) return std::clamp(yv.back(), 0.0, 1.0);
    const auto it = std::upper_bound(xq.begin(), xq.end(), q);
    const size_t i = static_cast<size_t>(std::distance(xq.begin(), it));
    const double t = (q - xq[i - 1]) / (xq[i] - xq[i - 1]);
    return std::clamp((1 - t) * yv[i - 1] + t * yv[i], 0.0, 1.0);
}

double AchillesCollision::FormFactor(PID pid, double q) const {
    const auto &ff = pid == PID::proton() ? m_ff_p : m_ff_n;
    if(ff.empty()) return 0;
    const double x = std::abs(q) / kFFdq;
    const size_t i = static_cast<size_t>(x);
    if(i + 1 >= ff.size()) return 0; // beyond 5 GeV: |F| < 1e-3 for any nucleus
    const double t = x - static_cast<double>(i);
    return (1 - t) * ff[i] + t * ff[i + 1];
}

double AchillesCollision::SigmaTot(Event &event, size_t i, size_t j) const {
    const double s = m_handler->TotalCrossSection(event, i, j);
    // a non-finite or negative value (e.g. the NASA parametrisation returns NaN when rounding
    // makes p_lab^2 < 0 at zero relative momentum) is treated as no interaction
    if(!std::isfinite(s) || s <= 0) return 0.0;
    return m_s.sigma_cap > 0 ? std::min(s, m_s.sigma_cap) : s;
}

double AchillesCollision::SigmaAtRest(Event &event, const Particle &proj) const {
    const auto key =
        std::make_tuple(static_cast<long>(proj.ID()), std::lround(proj.Momentum().P()), 0L);
    auto it = m_rest_cache.find(key);
    if(it != m_rest_cache.end()) return it->second;
    return m_rest_cache[key] = SigmaAtRestUncached(event, proj);
}

double AchillesCollision::SigmaAtRestUncached(Event &event, const Particle &proj_in) const {
    const Particle proj = proj_in; // copy: the hadron vector may reallocate below
    auto &hadrons = event.Hadrons();
    const size_t base = hadrons.size();
    hadrons.push_back(proj);
    hadrons.back().Status() = ParticleStatus::propagating;
    double sum = 0;
    for(PID pid : {PID::proton(), PID::neutron()}) {
        hadrons.emplace_back(pid, FourVector{ParticleInfo(pid).Mass(), 0, 0, 0}, proj.Position(),
                             ParticleStatus::background);
        sum += 0.5 * SigmaTot(event, base, base + 1);
        hadrons.pop_back();
    }
    hadrons.pop_back();
    return sum;
}

double AchillesCollision::SigmaBound(Event &event, const Particle &proj, double kmax) const {
    // cache per 1 MeV bin of |p|; the bound is evaluated at the upper edge of the bin and at
    // the bin centre and the larger is kept (sigma is smooth on that scale; the 10% safety
    // and the violation counter cover the rest)
    const double p = proj.Momentum().P();
    const auto key =
        std::make_tuple(static_cast<long>(proj.ID()), std::lround(p), std::lround(kmax));
    const double medium = m_medium ? m_s.medium_max : 1.0;
    auto it = m_bound_cache.find(key);
    if(it != m_bound_cache.end()) return medium * it->second;
    double best = 0;
    for(double dp : {-0.5, 0.0, 0.5}) {
        Particle q = proj;
        const double pn = std::max(1e-3, static_cast<double>(std::get<1>(key)) + dp);
        const ThreeVector dir = proj.Momentum().Vec3().Unit();
        const double m = proj.Info().Mass();
        q.SetMomentum(
            FourVector{std::sqrt(m * m + pn * pn), dir.Px() * pn, dir.Py() * pn, dir.Pz() * pn});
        best = std::max(best, SigmaBoundUncached(event, q, kmax));
    }
    m_bound_cache[key] = best;
    return medium * best;
}

double AchillesCollision::SigmaBoundUncached(Event &event, const Particle &proj_in,
                                             double kmax) const {
    const Particle proj = proj_in; // copy: the hadron vector may reallocate below
    // Scan target momenta |k| <= kmax over a set of directions; 10% safety. A miss is caught
    // by the violation counter.
    double best = 0;
    const std::array<PID, 2> pids{PID::proton(), PID::neutron()};
    auto &hadrons = event.Hadrons();
    const size_t base = hadrons.size();
    hadrons.push_back(proj);
    hadrons.back().Status() = ParticleStatus::propagating;
    const ThreeVector dir = proj.Momentum().Vec3().Unit();
    const ThreeVector perp = std::abs(dir.Px()) < 0.9 ? dir.Cross(ThreeVector{1, 0, 0}).Unit()
                                                      : dir.Cross(ThreeVector{0, 1, 0}).Unit();
    const size_t nk = kmax > 0 ? 12 : 1;
    for(PID pid : pids) {
        const double m = ParticleInfo(pid).Mass();
        for(size_t ik = 0; ik < nk; ++ik) {
            const double k =
                nk == 1 ? 0.0 : kmax * static_cast<double>(ik) / static_cast<double>(nk - 1);
            for(int c = -6; c <= 6; ++c) {
                const double cth = static_cast<double>(c) / 6.0;
                const double sth = std::sqrt(std::max(0.0, 1 - cth * cth));
                ThreeVector kv = dir * (k * cth) + perp * (k * sth);
                hadrons.emplace_back(
                    pid, FourVector{std::sqrt(m * m + k * k), kv.Px(), kv.Py(), kv.Pz()},
                    proj.Position(), ParticleStatus::background);
                best = std::max(best, SigmaTot(event, base, base + 1));
                hadrons.pop_back();
                if(k == 0) break;
            }
        }
    }
    hadrons.pop_back();
    return 1.10 * best;
}

void AchillesCollision::Propose(Event &event, Trial &trial, Context &ctx) {
    auto modes = m_handler->CrossSection(event, trial.proj, trial.target);
    double total = 0;
    for(const auto &m : modes) total += m.cross_section;
    if(modes.empty() || !std::isfinite(total) || total <= 0) {
        trial.factors.push_back({"channel", 0.0, {}});
        return;
    }
    auto mode = m_handler->SelectChannel(modes, ctx.uniform());
    trial.out_ids = mode.particles;
    const auto &p1 = event.Hadrons()[trial.proj];
    const auto &p2 = event.Hadrons()[trial.target];
    trial.outgoing = m_handler->GenerateMomentum(p1, p2, mode.particles, Random::Instance());
    bool ok = !trial.outgoing.empty();
    for(const auto &o : trial.outgoing)
        if(std::isnan(o.Momentum().Px())) ok = false;
    if(!ok) {
        trial.factors.push_back({"kinematics", 0.0, {}});
        return;
    }
    if(m_s.pauli) {
        double a = 1;
        for(const auto &o : trial.outgoing) {
            if(!o.Info().IsNucleon()) continue;
            a *= 1 - m_occ->f(o.Position().Magnitude(), o.ID(), o.Momentum().P());
        }
        trial.factors.push_back({"pauli", a, {}});
    }
    if(m_s.coherent) {
        // elastic: two nucleons out with the same species as in
        const bool elastic = trial.outgoing.size() == 2 && trial.outgoing[0].Info().IsNucleon() &&
                             trial.outgoing[1].Info().IsNucleon() &&
                             ((trial.outgoing[0].ID() == p1.ID() && trial.outgoing[1].ID() == p2.ID()) ||
                              (trial.outgoing[0].ID() == p2.ID() && trial.outgoing[1].ID() == p1.ID()));
        if(elastic) {
            // the projectile-like outgoing nucleon: the projectile's species, and the leading
            // one when both are (pp, nn). q = momentum transferred to the nucleus.
            const Particle *lead = nullptr;
            for(const auto &o : trial.outgoing) {
                if(o.ID() != p1.ID()) continue;
                if(!lead || o.Momentum().P() > lead->Momentum().P()) lead = &o;
            }
            const double q = (p1.Momentum().Vec3() - lead->Momentum().Vec3()).Magnitude();
            trial.factors.push_back({"coherent", 1 - CoherentVetoProb(p2.ID(), q), {}});
        }
    }
}

std::vector<size_t> AchillesCollision::Commit(Event &event, Trial &trial, Context &) {
    auto &parts = event.Hadrons();
    Particles initial{parts[trial.proj], parts[trial.target]};
    parts[trial.proj].Status() = ParticleStatus::interacted;
    parts[trial.target].Status() = ParticleStatus::interacted;
    initial[0].Status() = ParticleStatus::interacted;
    initial[1].Status() = ParticleStatus::interacted;
    const FourVector pin = parts[trial.proj].Momentum();
    Particles final_part;
    std::vector<size_t> added;
    for(auto out : trial.outgoing) {
        out.Status() = ParticleStatus::propagating;
        if(out.Info().IsNucleon()) out.SetFormationZone(pin, out.Momentum());
        parts.push_back(out);
        final_part.push_back(out);
        added.push_back(parts.size() - 1);
    }
    const auto pos = (initial[0].Position() + initial[1].Position()) / 2.0;
    event.History().AddVertex(pos, initial, final_part, EventHistory::StatusCode::cascade);
    return added;
}

// ==========================================================================================
// Factory
// ==========================================================================================
std::shared_ptr<Engine> MakeEngine(const YAML::Node &node, char track, Uniform01 uniform) {
    auto get = [&](const char *key, auto def) {
        using T = decltype(def);
        return node && node[key] ? node[key].as<T>() : def;
    };
    // Pauli: Step | FermiDirac | None (None: target Fermi motion kept, blocking off)
    const auto pauli = get("Pauli", std::string("Step"));
    if(pauli != "Step" && pauli != "FermiDirac" && pauli != "None")
        throw std::runtime_error("Sudakov: Pauli must be Step, FermiDirac or None");
    auto shape = pauli == "FermiDirac" ? Occupancy::Shape::FermiDirac : Occupancy::Shape::Step;
    auto occ = std::make_shared<Occupancy>(shape, get("PauliDelta", 0.05));
    AchillesCollisionSettings cs;
    cs.sigma_cap = get("SigmaCap", -1.0);
    cs.pauli = pauli != "None";
    cs.coherent = get("CoherentVeto", false);
    if(node && node["CoherentTable"]) {
        cs.coherent = true;
        cs.veto_q = node["CoherentTable"]["q"].as<std::vector<double>>();
        cs.veto_v = node["CoherentTable"]["v"].as<std::vector<double>>();
        if(cs.veto_q.size() != cs.veto_v.size() || cs.veto_q.size() < 2 ||
           !std::is_sorted(cs.veto_q.begin(), cs.veto_q.end()))
            throw std::runtime_error("Sudakov: CoherentTable needs sorted q and v of equal length");
    }
    cs.medium_max = get("MediumMax", 1.0);
    auto coll = std::make_shared<AchillesCollision>(occ, cs);
    EngineSettings es;
    es.max_trials_per_event = get("MaxTrials", size_t{2000000});
    // one attempt per pair (Achilles' prescription and the classical-profile rule of Step 1-2
    // section 5); Retry re-tests the pair at each later majorant point
    const auto attempt = get("PauliAttempt", std::string("Single"));
    if(attempt != "Single" && attempt != "Retry")
        throw std::runtime_error("Sudakov: PauliAttempt must be Single or Retry");
    es.single_attempt = attempt == "Single";
    std::vector<double> scales;
    if(node && node["Variations"]) {
        for(const auto &v : node["Variations"]) {
            es.variations.push_back(v["Name"].as<std::string>());
            scales.push_back(v["SigmaScale"].as<double>());
        }
    }
    std::shared_ptr<RateModel> rate;
    if(track == 'A') {
        DensityRateSettings s;
        s.shells = get("Shells", size_t{16});
        s.rmax = get("RMax", -1.0);
        s.headroom = get("Headroom", 1.0);
        s.fold_width = get("FoldWidth", 0.0);
        s.finite_A = get("FiniteA", false);
        s.sigma_scale_vars = scales;
        if(node && node["Hole"]) {
            s.hole_d = node["Hole"]["d"].as<std::vector<double>>();
            s.hole_g = node["Hole"]["g"].as<std::vector<double>>();
            s.hole_gmax = *std::max_element(s.hole_g.begin(), s.hole_g.end());
            s.hole_gmax = std::max(1.0, s.hole_gmax);
        }
        rate = std::make_shared<DensityRate>(s, coll, occ);
    } else {
        ConfigurationRateSettings s;
        const auto law = get("Law", std::string("Eikonal"));
        s.law = law == "Kappa1" ? ConfigurationRateSettings::Law::Kappa1
                                : ConfigurationRateSettings::Law::Eikonal;
        const auto prof = get("Profile", std::string("Gauss"));
        s.profile = prof == "Atom"  ? ConfigurationRateSettings::Profile::Atom
                    : prof == "Box" ? ConfigurationRateSettings::Profile::Box
                                    : ConfigurationRateSettings::Profile::Gauss;
        s.x0 = get("X0", 2.6);
        s.headroom = get("Headroom", 1.0);
        s.sigma_scale_vars = scales;
        rate = std::make_shared<ConfigurationRate>(s, coll);
        // Track B uses the configuration nucleons; the occupancy still needs densities
    }
    return std::make_shared<Engine>(rate, coll, es, std::move(uniform));
}

} // namespace achilles::sudakov
