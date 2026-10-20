// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/ParticleInfo.hh"
#include "Achilles/NuclearMass.hh"
#include "Achilles/System.hh"

#include <cstdlib>
#include <iostream>
#include <memory>

using achilles::ParticleInfo;
using achilles::ParticleInfoEntry;
using achilles::PID;

// achilles::Database::ParticleDB achilles::Database::particleDB{};

// PID::PID(const long int &id_, bool validate) : id(PID::undefined()) {
//     if(!validate)
//          id = std::abs(id_);
//     else
//          valid(std::abs(id_));
//
//     id = std::abs(id_);
//     // Check if fundamental SM particles
//     if(id > 0 && id < 30) {
//     }
// }

bool PID::valid_nucleus() const {
    // Nuclear codes are +/- 10LZZZAAAI, so the leading two digits must be 10.
    if(Abs_() / 1000000000 != 1) return false;
    return NuclearA() >= 1 && NuclearL() + NuclearZ() <= NuclearA();
}

std::shared_ptr<ParticleInfoEntry> ParticleInfo::Lookup(PID id) {
    auto it = particleDB.find(id);
    if(it != particleDB.end()) return it->second;
    if(id.valid_nucleus()) return BuildNucleusEntry(id);
    throw std::runtime_error(fmt::format("Invalid PID: id={}", id.AsInt()));
}

std::shared_ptr<ParticleInfoEntry> ParticleInfo::BuildNucleusEntry(PID id) {
    const int L = id.NuclearL(), Z = id.NuclearZ(), A = id.NuclearA(), I = id.NuclearI();
    const std::string name = NuclearName(Z, A, L, I);

    auto entry = std::make_shared<ParticleInfoEntry>(id, NuclearMass(Z, A), 0.0, 3 * Z, 0,
                                                     std::abs(A - 2 * Z), A % 2, 1, 0, true, true,
                                                     name, "anti-" + name);
    particleDB.emplace(id, entry);
    nameToPID.emplace(name, id);
    return entry;
}

void achilles::ParticleInfo::BuildDatabase(const std::string &datafile) {
    YAML::Node particleYAML = YAML::LoadFile(Filesystem::FindFile(datafile, "ParticleInfo"));
    auto particles = particleYAML["Particles"];
    for(auto particle : particles) {
        auto entry =
            std::make_shared<ParticleInfoEntry>(particle["Particle"].as<ParticleInfoEntry>());
        particleDB.emplace(entry->id, entry);
        nameToPID.emplace(entry->idname, entry->id);
    }
    PrintDatabase();
}

void achilles::ParticleInfo::PrintDatabase() {
    fmt::print("{:>10s} {:<20s} {:<20s} {:^6s} {:^10s}    {:^10s}\n", "PID", "Name", "Anti-name",
               "Stable", "Mass (MeV)", "Width (MeV)");
    for(const auto &part : particleDB) { std::cout << *(part.second) << "\n"; }
}

ParticleInfoEntry::ParticleInfoEntry(PID id_, double mass_, double width_, int icharge_,
                                     int strong_, int isospin_, int spin_, int stable_,
                                     int majorana_, bool massive_, bool hadron_,
                                     std::string idname_, std::string antiname_)
    : id(id_), mass(mass_), hmass(mass_), width(width_), icharge(icharge_), strong(strong_),
      isospin(isospin_), spin(spin_), stable(stable_), majorana(majorana_), massive(massive_),
      hadron(hadron_), idname(std::move(idname_)), antiname(std::move(antiname_)) {}

std::ostream &operator<<(std::ostream &os, const ParticleInfoEntry &entry) {
    os << entry.ToString();
    return os;
}

ParticleInfo::ParticleDB ParticleInfo::particleDB;
std::map<std::string, achilles::PID> ParticleInfo::nameToPID;

bool ParticleInfo::IsNucleon() const noexcept {
    return info->id == PID::proton() || info->id == PID::neutron();
}

bool ParticleInfo::IsDelta() const noexcept {
    return info->id == PID::deltapp() || info->id == PID::deltap() || info->id == PID::delta0() ||
           info->id == PID::deltam();
}

bool ParticleInfo::IsPion() const noexcept {
    return info->id == PID::pionp() || info->id == PID::pion0() || info->id == -PID::pionp();
}

bool ParticleInfo::IsBaryon() const noexcept {
    if(IntID() % 10000 < 1000) return false;
    return true;
}

bool ParticleInfo::IsBHadron() const noexcept {
    if(IntID() < 100) return false;
    if(IntID() - 100 * IntID() / 100 < 10) return false;
    if((IntID() - 100 * IntID() / 100) / 10 == 5) return true;
    if((IntID() - 1000 * IntID() / 1000) / 100 == 5) return true;
    if((IntID() - 10000 * IntID() / 10000) / 1000 == 5) return true;
    return false;
}

bool ParticleInfo::IsCHadron() const noexcept {
    if(IntID() < 100) return false;
    if(IntID() - 100 * IntID() / 100 < 10) return false;
    if((IntID() - 100 * IntID() / 100) / 10 == 4) return true;
    if((IntID() - 1000 * IntID() / 1000) / 100 == 4) return true;
    if((IntID() - 10000 * IntID() / 10000) / 1000 == 4) return true;
    return false;
}

bool ParticleInfo::IsSHadron() const noexcept {
    if(IntID() < 100) return false;
    if(IntID() - 100 * IntID() / 100 < 10) return false;
    if((IntID() - 100 * IntID() / 100) / 10 == 3) return true;
    if((IntID() - 1000 * IntID() / 1000) / 100 == 3) return true;
    if((IntID() - 10000 * IntID() / 10000) / 1000 == 3) return true;
    return false;
}

size_t ParticleInfo::NSpins() const noexcept {
    if(IsFermion()) {
        if(IntSpin() == 1)
            return 2;
        else
            return 4;
    } else {
        if(IsVector())
            return IsMassive() ? 3 : 2;
        else
            return 1;
    }
}

double ParticleInfo::GenerateLifeTime() const {
    throw std::runtime_error("Not Implemented Yet");
    return 0.0;
}

bool ParticleInfo::IsStable() const noexcept {
    if(info->stable == 0) return false;
    if(info->stable == 1) return true;
    if(info->stable == 2 && !IsAnti()) return true;
    if(info->stable == 3 && IsAnti()) return true;
    return false;
}
