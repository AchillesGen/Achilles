// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

// achilles-cascade-study <runcard> <output.csv> [seed]
//
// Event-by-event diagnostics of the cascade for the validation study of the curved-path
// algorithms. Uses the same run card as achilles-cascade (Nucleus, Cascade, PID) plus
//   Study: {NEvents: <events per momentum>, Momenta: [...], Survival: [times], NPauli: n}
// Every event is written (hit or not) with the first-collision record, the leading final
// nucleon, the algorithm counters and the CPU time. Transparency events also carry the
// struck nucleon's initial (Fermi) momentum for missing-momentum cuts.

#include "Achilles/Cascade.hh"
#include "Achilles/CurvedPath.hh"
#include "Achilles/Event.hh"
#include "Achilles/EventHistory.hh"
#include "Achilles/Logging.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Potential.hh"
#include "Achilles/Random.hh"
#include "Achilles/RunCascade.hh"
#include "Achilles/Settings.hh"

#include "fmt/format.h"
#include "fmt/os.h"
#include "spdlog/spdlog.h"

#include <chrono>
#include <cmath>
#include <string>

using namespace achilles;

namespace {

struct Kick {
    ThreeVector q0, p_initial, p_kick;
    size_t idx{};
};

Kick KickInternal(Event &event, PID pid, double mom, std::shared_ptr<Nucleus> nuc,
                  bool formation_zone) {
    // As CascadeTest::InitTransparency (external: false), recording the initial momentum
    const double costheta = Random::Instance().Uniform(-1.0, 1.0);
    const double sintheta = std::sqrt(1 - costheta * costheta);
    const double phi = Random::Instance().Uniform(0.0, 2 * M_PI);
    std::vector<FourVector> momdummy = {};
    event = Event(nuc, momdummy, 1.);
    std::vector<size_t> candidates;
    for(size_t i = 0; i < event.Hadrons().size(); ++i)
        if(event.Hadrons()[i].ID() == pid) candidates.push_back(i);
    const size_t idx = candidates[Random::Instance().Uniform(0ul, candidates.size() - 1)];
    auto &kicked = event.Hadrons()[idx];
    const double mass = kicked.Info().Mass();
    FourVector kick{std::sqrt(mom * mom + mass * mass), mom * sintheta * std::cos(phi),
                    mom * sintheta * std::sin(phi), mom * costheta};
    Kick record{kicked.Position(), kicked.Momentum().Vec3(), kick.Vec3(), idx};
    if(formation_zone) kicked.SetFormationZone(kicked.Momentum(), kick);
    kicked.Status() = ParticleStatus::internal_test;
    kicked.SetMomentum(kick);
    return record;
}

} // namespace

int main(int argc, char *argv[]) {
    if(argc < 3) {
        fmt::print("Usage: achilles-cascade-study <runcard> <output.csv> [seed]\n");
        return 1;
    }
    CreateLogger(3, 3, 5);
    Settings setting{argv[1], "data/hA_Rules.yml"};
    unsigned int seed = 123456789;
    if(setting["Initialize"]["seed"]) seed = setting["Initialize"]["seed"].as<unsigned int>();
    if(argc > 3) seed = static_cast<unsigned int>(std::stoul(argv[3]));
    Random::Instance().Seed(seed);

    const auto mode = setting.GetAs<CascadeMode>("Cascade/Mode");
    auto cascade = setting.GetAs<Cascade>("Cascade");
    const auto params = setting.GetAs<std::map<std::string, double>>("Cascade/Params");
    const auto pid = setting.GetAs<PID>("PID");
    auto nuc = std::make_shared<Nucleus>(setting.GetAs<Nucleus>("Nucleus"));
    auto potential = PotentialFactory::Initialize(
        setting.GetAs<std::string>("Nucleus/Potential/Name"), nuc, setting["Nucleus/Potential"]);
    const auto study_node = setting["Study"];
    if(study_node["TabulatePotential"] && study_node["TabulatePotential"].as<bool>()) {
        // all algorithms propagate in the same bicubic table of the nuclear potential
        nuc->SetPotential(std::make_unique<TabulatedPotential>(*potential));
    } else {
        nuc->SetPotential(std::move(potential));
    }

    const auto study = setting["Study"];
    const auto nevents = study["NEvents"].as<size_t>();
    const auto momenta = study["Momenta"].as<std::vector<double>>();
    std::vector<double> survival_times;
    if(study["Survival"]) survival_times = study["Survival"].as<std::vector<double>>();
    const size_t npauli = study["NPauli"] ? study["NPauli"].as<size_t>() : 0;
    const bool formation_zone = study["FormationZone"] ? study["FormationZone"].as<bool>() : true;
    // Diagnostic: evolve copies of each event Repeat times and record the surviving fraction
    // at the Survival times (columns R0..), to compare with SurvivalProfile event by event
    const size_t repeat = study["Repeat"] ? study["Repeat"].as<size_t>() : 0;

    auto out = fmt::output_file(argv[2]);
    out.print("mom,hit,ncoll,first,t1,sep1,b1,z1,captured,lead_px,lead_py,lead_pz,lead_E,lead_id,"
              "kick_px,kick_py,kick_pz,pi_x,pi_y,pi_z,q0x,q0y,q0z,ncand,nacc,nblock,nviol,"
              "maxratio,weight,nsub,nmin,npre,emis,usec,h0,pref_x,pref_y,pref_z,nunres,ncontact");
    for(size_t k = 0; k < survival_times.size(); ++k) out.print(",S{}", k);
    if(repeat > 0)
        for(size_t k = 0; k < survival_times.size(); ++k) out.print(",R{}", k);
    out.print("\n");
    fmt::print("radius = {} fm, algorithm run, seed {}\n", nuc->Radius(), seed);

    size_t nfailed = 0;
    for(const double mom : momenta) {
        for(size_t iev = 0; iev < nevents; ++iev) {
            Event event;
            Kick kick;
            if(mode == CascadeMode::Transparency) {
                kick = KickInternal(event, pid, mom, nuc, formation_zone);
            } else {
                CascadeTest::InitCrossSection(event, pid, mom, params.at("radius"), nuc);
                kick.idx = event.Hadrons().size() - 1;
                kick.q0 = event.Hadrons().back().Position();
                kick.p_kick = event.Hadrons().back().Momentum().Vec3();
            }
            cascade.SetKicked(kick.idx);
            // Reference without collisions: conserved energy and asymptotic momentum of the
            // kicked particle on its own path (mean field on: curved; off: the kick itself)
            double h0 = 0;
            ThreeVector pref = kick.p_kick;
            {
                const auto &kp = event.Hadrons()[kick.idx];
                const double mass = kp.Info().Mass();
                const bool pot = cascade.UsePotentialProp() && nuc->GetPotential() &&
                                 kp.Info().IsNucleon() && mode == CascadeMode::Transparency;
                CurvedPath path(0, kp.Position(), kp.Momentum().Vec3(), mass,
                                pot ? nuc->GetPotential() : nullptr, 0.04, nuc->Radius(), false);
                h0 = path.Energy();
                if(pot && h0 > mass && !path.Trapped()) {
                    const double tex = path.ExitTime();
                    if(std::isfinite(tex)) {
                        const auto p = path.Momentum(tex);
                        pref = std::sqrt(h0 * h0 - mass * mass) * p / p.P();
                    }
                }
            }
            std::vector<double> surv;
            if(!survival_times.empty()) {
                // Status must be propagating for SurvivalProfile's track; it reads positions only
                surv = cascade.SurvivalProfile(event, nuc.get(), kick.idx, survival_times, npauli);
            }
            std::vector<double> rep(survival_times.size(), 0);
            for(size_t r = 0; r < repeat; ++r) {
                Event copy = event;
                cascade.SetKicked(kick.idx);
                try {
                    cascade.Evolve(copy, nuc.get());
                } catch(const std::exception &) { cascade.Reset(); }
                const auto &dg = cascade.GetDiagnostics();
                for(size_t k = 0; k < survival_times.size(); ++k)
                    rep[k] += (!dg.first_set || dg.first_time > survival_times[k]) / static_cast<double>(repeat);
                cascade.Reset();
                cascade.SetKicked(kick.idx);
            }
            const auto start = std::chrono::steady_clock::now();
            bool failed = false;
            try {
                cascade.Evolve(event, nuc.get());
            } catch(const std::exception &e) {
                // e.g. Validate: a particle still propagating after the maximal number of steps
                failed = true;
                cascade.Reset();
                ++nfailed;
            }
            const auto stop = std::chrono::steady_clock::now();
            const double usec =
                std::chrono::duration<double, std::micro>(stop - start).count();
            const auto &diag = cascade.GetDiagnostics();

            // Leading final-state particle of the kicked species (largest momentum)
            const Particle *lead = nullptr;
            size_t ncaptured = 0;
            for(const auto &part : event.Hadrons()) {
                if(part.Status() == ParticleStatus::captured) ++ncaptured;
                if(part.Status() != ParticleStatus::final_state || part.ID() != pid) continue;
                if(!lead || part.Momentum().P() > lead->Momentum().P()) lead = &part;
            }
            const FourVector lp = lead ? lead->Momentum() : FourVector{};
            out.print("{},{},{},{},{:.9g},{:.9g},{:.9g},{:.9g},{},{:.9g},{:.9g},{:.9g},{:.9g},{},"
                      "{:.9g},{:.9g},{:.9g},{:.9g},{:.9g},{:.9g},{:.9g},{:.9g},{:.9g},{},{},{},{},"
                      "{:.6g},{:.9g},{},{},{},{:.6g},{:.4g}",
                      mom, failed ? -1 : (event.History().size() > 0 ? 1 : 0), diag.collisions,
                      diag.first_set ? 1 : 0, diag.first_time, diag.first_separation, diag.first_b,
                      diag.first_z, ncaptured, lp.Px(), lp.Py(), lp.Pz(), lp.E(),
                      lead ? lead->ID().AsInt() : 0L, kick.p_kick.Px(), kick.p_kick.Py(),
                      kick.p_kick.Pz(), kick.p_initial.Px(), kick.p_initial.Py(),
                      kick.p_initial.Pz(), kick.q0.Px(), kick.q0.Py(), kick.q0.Pz(),
                      diag.candidates, diag.accepted, diag.pauli_blocked, diag.violations,
                      diag.max_ratio, diag.weight, diag.subdivisions, diag.minima,
                      diag.prescreened, diag.energy_mismatch, usec);
            out.print(",{:.9g},{:.9g},{:.9g},{:.9g},{},{}", h0, pref.Px(), pref.Py(), pref.Pz(),
                      diag.unresolved, diag.contacts);
            for(const double s : surv) out.print(",{:.9g}", s);
            if(repeat > 0)
                for(const double s : rep) out.print(",{:.9g}", s);
            out.print("\n");
            cascade.Reset();
        }
        fmt::print("p = {} MeV done, failed events so far: {}\n", mom, nfailed);
    }
    return 0;
}
