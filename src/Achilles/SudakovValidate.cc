// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later
//
// achilles-sudakov-validate: drives the Sudakov engine in forward (no-deflection) mode and
// writes per-event CSV output. The oracles (closed forms, enumeration, KS/chi2) live in
// scripts/sudakov_validation/analyse.py, so the code under test and the reference are separate.
//
// usage: achilles-sudakov-validate <test> [key=value ...] > out.csv
//   tests: sphere-ext, sphere-int, toy, qmc-ext, qmc-int
//   keys : track=A|B sigma=40 headroom=1 nev=100000 seed=1 law=Kappa1|Eikonal
//          profile=Atom|Box|Gauss x0=2.6 vars=0.8,1.2 rhohat=-1 finiteA=0 fold=0
//          configs=file density=file decorrelate=0 pauli=1 toy=2|3|4 internal=0 pmom=1000

#include "Achilles/Event.hh"
#include "Achilles/Particle.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/Sudakov.hh"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace achilles;
using namespace achilles::sudakov;

namespace {
std::map<std::string, std::string> g_args;
std::string Arg(const std::string &k, const std::string &def) {
    auto it = g_args.find(k);
    return it == g_args.end() ? def : it->second;
}
double ArgD(const std::string &k, double def) {
    auto it = g_args.find(k);
    return it == g_args.end() ? def : std::stod(it->second);
}
std::vector<double> ArgV(const std::string &k) {
    std::vector<double> out;
    std::stringstream ss(Arg(k, ""));
    std::string tok;
    while(std::getline(ss, tok, ','))
        if(!tok.empty()) out.push_back(std::stod(tok));
    return out;
}

std::mt19937_64 g_rng;
FILE *g_out = stdout;
double U() {
    // (0, 1)
    return (static_cast<double>(g_rng() >> 11) + 0.5) * 0x1.0p-53;
}

constexpr double A12 = 12, RMS = 2.33;
const double Rsph = std::sqrt(5.0 / 3.0) * RMS;
const double rho0 = A12 / (4.0 / 3.0 * M_PI * Rsph * Rsph * Rsph);

struct Config {
    std::vector<ThreeVector> pos;
    std::vector<int> iso;
    double w;
};
std::vector<Config> g_configs;
std::vector<double> g_cdf; // weights CDF
void LoadConfigs(const std::string &fn) {
    std::ifstream in(fn);
    if(!in) throw std::runtime_error("cannot open configs " + fn);
    size_t N, A;
    in >> N >> A;
    g_configs.resize(N);
    for(auto &c : g_configs) {
        c.pos.resize(A);
        c.iso.resize(A);
        for(size_t j = 0; j < A; ++j) {
            double x, y, z;
            in >> c.iso[j] >> x >> y >> z;
            c.pos[j] = {x, y, z};
        }
        in >> c.w;
    }
    g_cdf.resize(N);
    double acc = 0;
    for(size_t i = 0; i < N; ++i) g_cdf[i] = (acc += g_configs[i].w);
    for(auto &c : g_cdf) c /= acc;
}
size_t PickConfig() {
    return static_cast<size_t>(std::lower_bound(g_cdf.begin(), g_cdf.end(), U()) - g_cdf.begin());
}
ThreeVector RandomRotate(const ThreeVector &v, const std::array<double, 9> &R) {
    return {R[0] * v.Px() + R[1] * v.Py() + R[2] * v.Pz(),
            R[3] * v.Px() + R[4] * v.Py() + R[5] * v.Pz(),
            R[6] * v.Px() + R[7] * v.Py() + R[8] * v.Pz()};
}
std::array<double, 9> RandomRotation() {
    // uniform quaternion
    double q[4];
    double n = 0;
    std::normal_distribution<double> g;
    for(auto &x : q) {
        x = g(g_rng);
        n += x * x;
    }
    n = std::sqrt(n);
    for(auto &x : q) x /= n;
    const double a = q[0], b = q[1], c = q[2], d = q[3];
    return {a * a + b * b - c * c - d * d, 2 * (b * c - a * d),
            2 * (b * d + a * c),           2 * (b * c + a * d),
            a * a - b * b + c * c - d * d, 2 * (c * d - a * b),
            2 * (b * d - a * c),           2 * (c * d + a * b),
            a * a - b * b - c * c + d * d};
}

// radial density table (r, rho_p, rho_n)
std::vector<double> g_tr, g_tp, g_tn;
void LoadDensity(const std::string &fn) {
    std::ifstream in(fn);
    if(!in) throw std::runtime_error("cannot open density " + fn);
    double r, p, n;
    while(in >> r >> p >> n) {
        g_tr.push_back(r);
        g_tp.push_back(p);
        g_tn.push_back(n);
    }
}
double Interp(const std::vector<double> &y, double r) {
    if(r >= g_tr.back()) return 0;
    auto it = std::upper_bound(g_tr.begin(), g_tr.end(), r);
    const size_t i = static_cast<size_t>(it - g_tr.begin());
    if(i == 0) return y[0];
    const double f = (r - g_tr[i - 1]) / (g_tr[i] - g_tr[i - 1]);
    return (1 - f) * y[i - 1] + f * y[i];
}

Particle Projectile(const ThreeVector &pos, ParticleStatus status) {
    const double p = ArgD("pmom", 1000.0);
    const double m = ParticleInfo(PID::proton()).Mass();
    return Particle(PID::proton(), FourVector{std::sqrt(m * m + p * p), 0, 0, p}, pos, status);
}

std::shared_ptr<Engine> MakeForwardEngine(char track) {
    const double sigma = ArgD("sigma", 40);
    const auto vars = ArgV("vars");
    const double pa = ArgD("pauli", 1.0);
    auto coll = std::make_shared<ForwardCollision>(sigma, pa, ArgV("pvars"));
    EngineSettings es;
    for(size_t i = 0; i < std::max(vars.size(), ArgV("pvars").size()); ++i)
        es.variations.push_back("v" + std::to_string(i));
    std::shared_ptr<RateModel> rate;
    if(track == 'A') {
        DensityRateSettings s;
        s.headroom = ArgD("headroom", 1.0);
        s.shells = static_cast<size_t>(ArgD("shells", 16));
        s.sigma_scale_vars = vars;
        s.rho_hat_override = ArgD("rhohat", -1);
        s.finite_A = ArgD("finiteA", 0) > 0;
        s.fold_width = ArgD("fold", 0);
        auto occ = std::make_shared<Occupancy>(Occupancy::Shape::Step, 0.05);
        // targets at rest: zero occupancy densities
        occ->SetDensities([](double) { return 0.0; }, [](double) { return 0.0; });
        auto dr = std::make_shared<DensityRate>(s, coll, occ);
        if(!g_tr.empty()) {
            dr->SetDensities([](double r) { return Interp(g_tp, r); },
                             [](double r) { return Interp(g_tn, r); }, A12, g_tr.back());
        } else {
            dr->SetDensities([](double r) { return r < Rsph ? rho0 / 2 : 0.0; },
                             [](double r) { return r < Rsph ? rho0 / 2 : 0.0; }, A12, Rsph);
        }
        rate = dr;
    } else {
        ConfigurationRateSettings s;
        s.law = Arg("law", "Kappa1") == "Kappa1" ? ConfigurationRateSettings::Law::Kappa1
                                                 : ConfigurationRateSettings::Law::Eikonal;
        const auto pr = Arg("profile", "Box");
        s.profile = pr == "Atom"  ? ConfigurationRateSettings::Profile::Atom
                    : pr == "Box" ? ConfigurationRateSettings::Profile::Box
                                  : ConfigurationRateSettings::Profile::Gauss;
        s.x0 = ArgD("x0", 2.6);
        s.headroom = ArgD("headroom", 1.0);
        s.sigma_scale_vars = vars;
        rate = std::make_shared<ConfigurationRate>(s, coll);
    }
    return std::make_shared<Engine>(rate, coll, es, U);
}

void Header(const char *cols, const Engine &eng) {
    std::fprintf(g_out, "%s", cols);
    for(size_t v = 0; v < eng.Settings().variations.size(); ++v) std::fprintf(g_out, ",w%zu", v);
    std::fprintf(g_out, "\n");
}

void PrintWeights(const Engine &eng) {
    for(double w : eng.GetContext().weights) std::fprintf(g_out, ",%.10g", w);
}

void PrintDiag(const Engine &eng) {
    const auto &d = eng.GetDiagnostics();
    std::fprintf(stderr,
                 "events=%zu trials=%zu accepted=%zu violations=%zu max_ratio=%.4f capped=%zu\n",
                 d.events, d.trials, d.accepted, d.violations, d.max_ratio, d.capped);
}

// ----------------------------------------------------------------------------------------
void SphereExt(char track) {
    auto eng = MakeForwardEngine(track);
    const size_t nev = static_cast<size_t>(ArgD("nev", 100000));
    std::poisson_distribution<int> pois(A12);
    Header("b,nint,depth", *eng);
    for(size_t ev = 0; ev < nev; ++ev) {
        const double b = Rsph * std::sqrt(U()), ph = 2 * M_PI * U();
        const double zin = -std::sqrt(Rsph * Rsph - b * b);
        Event event;
        event.Hadrons().push_back(Projectile({b * std::cos(ph), b * std::sin(ph), -Rsph - 10.0},
                                             ParticleStatus::external_test));
        if(track == 'B') {
            const int n = pois(g_rng);
            for(int j = 0; j < n; ++j) {
                ThreeVector x;
                do {
                    x = {Rsph * (2 * U() - 1), Rsph * (2 * U() - 1), Rsph * (2 * U() - 1)};
                } while(x.Magnitude2() >= Rsph * Rsph);
                event.Hadrons().emplace_back(
                    j % 2 ? PID::proton() : PID::neutron(),
                    FourVector{ParticleInfo(PID::proton()).Mass(), 0, 0, 0}, x,
                    ParticleStatus::background);
            }
        }
        eng->Evolve(event, nullptr);
        const auto &cols = eng->GetContext().collisions;
        const double depth = cols.empty() ? -1 : cols.front().point.Pz() - zin;
        std::fprintf(g_out, "%.6f,%zu,%.6f", b, cols.size(), depth);
        PrintWeights(*eng);
        std::fprintf(g_out, "\n");
    }
    PrintDiag(*eng);
}

void SphereInt(char track) {
    auto eng = MakeForwardEngine(track);
    const size_t nev = static_cast<size_t>(ArgD("nev", 100000));
    Header("nint", *eng);
    for(size_t ev = 0; ev < nev; ++ev) {
        ThreeVector x;
        do {
            x = {Rsph * (2 * U() - 1), Rsph * (2 * U() - 1), Rsph * (2 * U() - 1)};
        } while(x.Magnitude2() >= Rsph * Rsph);
        Event event;
        event.Hadrons().push_back(Projectile(x, ParticleStatus::internal_test));
        eng->Evolve(event, nullptr);
        std::fprintf(g_out, "%zu", eng->GetContext().collisions.size());
        PrintWeights(*eng);
        std::fprintf(g_out, "\n");
    }
    PrintDiag(*eng);
}

void Toy() {
    auto eng = MakeForwardEngine('B');
    const int which = static_cast<int>(ArgD("toy", 3));
    const bool internal = ArgD("internal", 0) > 0;
    std::vector<ThreeVector> pts;
    if(which == 2) pts = {{0.3, 0, -0.4}, {0.5, 0.2, 0.4}};
    if(which == 3) pts = {{0.3, 0, -1.0}, {0.8, 0, 0}, {0.2, 0.1, 0.5}};
    if(which == 4) pts = {{0, 0.3, -1.5}, {0.6, 0, -0.3}, {0.1, -0.4, 0.3}, {0.9, 0.5, 1.2}};
    const size_t nev = static_cast<size_t>(ArgD("nev", 200000));
    std::fprintf(g_out, "mask,first,tfirst\n");
    for(size_t ev = 0; ev < nev; ++ev) {
        Event event;
        event.Hadrons().push_back(
            Projectile({0, 0, internal ? 0.0 : -20.0},
                       internal ? ParticleStatus::internal_test : ParticleStatus::external_test));
        for(const auto &p : pts)
            event.Hadrons().emplace_back(PID::neutron(),
                                         FourVector{ParticleInfo(PID::neutron()).Mass(), 0, 0, 0},
                                         p, ParticleStatus::background);
        eng->Evolve(event, nullptr);
        unsigned mask = 0;
        int first = -1;
        double tf = -1;
        for(const auto &c : eng->GetContext().collisions) {
            mask |= 1u << (c.target - 1);
            if(first < 0) {
                first = static_cast<int>(c.target) - 1;
                tf = c.point.Pz();
            }
        }
        std::fprintf(g_out, "%u,%d,%.6f\n", mask, first, tf);
    }
    PrintDiag(*eng);
}

void FillConfig(Event &event, bool decorrelate) {
    const auto R = RandomRotation();
    const size_t A = g_configs.front().pos.size();
    if(!decorrelate) {
        const auto &c = g_configs[PickConfig()];
        for(size_t j = 0; j < A; ++j)
            event.Hadrons().emplace_back(c.iso[j] > 0 ? PID::neutron() : PID::proton(),
                                         FourVector{ParticleInfo(PID::proton()).Mass(), 0, 0, 0},
                                         RandomRotate(c.pos[j], R), ParticleStatus::background);
    } else {
        for(size_t j = 0; j < A; ++j) {
            const auto &c = g_configs[PickConfig()];
            const size_t s = static_cast<size_t>(U() * static_cast<double>(A)) % A;
            event.Hadrons().emplace_back(c.iso[s] > 0 ? PID::neutron() : PID::proton(),
                                         FourVector{ParticleInfo(PID::proton()).Mass(), 0, 0, 0},
                                         RandomRotate(c.pos[s], R), ParticleStatus::background);
        }
    }
}

void QmcExt(char track) {
    auto eng = MakeForwardEngine(track);
    const size_t nev = static_cast<size_t>(ArgD("nev", 100000));
    const bool dec = ArgD("decorrelate", 0) > 0;
    const double bmax = ArgD("bmax", 8.0);
    Header("b,nint", *eng);
    for(size_t ev = 0; ev < nev; ++ev) {
        const double b = bmax * std::sqrt(U()), ph = 2 * M_PI * U();
        Event event;
        event.Hadrons().push_back(
            Projectile({b * std::cos(ph), b * std::sin(ph), -12.0}, ParticleStatus::external_test));
        if(track == 'B') FillConfig(event, dec);
        eng->Evolve(event, nullptr);
        std::fprintf(g_out, "%.5f,%zu", b, eng->GetContext().collisions.size());
        PrintWeights(*eng);
        std::fprintf(g_out, "\n");
    }
    PrintDiag(*eng);
}

void QmcInt(char track) {
    auto eng = MakeForwardEngine(track);
    const size_t nev = static_cast<size_t>(ArgD("nev", 100000));
    const bool dec = ArgD("decorrelate", 0) > 0;
    Header("nint", *eng);
    for(size_t ev = 0; ev < nev; ++ev) {
        Event event;
        // knockout nucleon from a configuration (Track B) or from the density (Track A)
        Event cfg;
        FillConfig(cfg, dec);
        const size_t j0 = static_cast<size_t>(U() * static_cast<double>(cfg.Hadrons().size())) %
                          cfg.Hadrons().size();
        event.Hadrons().push_back(
            Projectile(cfg.Hadrons()[j0].Position(), ParticleStatus::internal_test));
        if(track == 'B')
            for(size_t j = 0; j < cfg.Hadrons().size(); ++j)
                if(j != j0) event.Hadrons().push_back(cfg.Hadrons()[j]);
        eng->Evolve(event, nullptr);
        std::fprintf(g_out, "%zu", eng->GetContext().collisions.size());
        PrintWeights(*eng);
        std::fprintf(g_out, "\n");
    }
    PrintDiag(*eng);
}
} // namespace

int main(int argc, char **argv) {
    if(argc < 2) {
        std::cerr << "usage: achilles-sudakov-validate <test> [key=value ...]\n";
        return 1;
    }
    const std::string test = argv[1];
    for(int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        auto eq = a.find('=');
        if(eq != std::string::npos) g_args[a.substr(0, eq)] = a.substr(eq + 1);
    }
    g_rng.seed(static_cast<unsigned long>(ArgD("seed", 1)));
    if(!Arg("out", "").empty()) g_out = std::fopen(Arg("out", "").c_str(), "w");
    if(!Arg("configs", "").empty()) LoadConfigs(Arg("configs", ""));
    if(!Arg("density", "").empty()) LoadDensity(Arg("density", ""));
    const char track = Arg("track", "A")[0];
    if(test == "sphere-ext")
        SphereExt(track);
    else if(test == "sphere-int")
        SphereInt(track);
    else if(test == "toy")
        Toy();
    else if(test == "qmc-ext")
        QmcExt(track);
    else if(test == "qmc-int")
        QmcInt(track);
    else {
        std::cerr << "unknown test " << test << "\n";
        return 1;
    }
    if(g_out != stdout) std::fclose(g_out);
    return 0;
}
