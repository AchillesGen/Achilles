// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Achilles/ParticleInfo.hh"
#include "Achilles/System.hh"

#include <algorithm>
#include <cmath>
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
    long int L = id / 10000000 - 100;
    long int Z = (id % 10000000) / 10000;
    long int A = (id % 10000) / 10;
    return L + Z < A;
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
    if(particleYAML["Containers"]) {
        for(auto container : particleYAML["Containers"]) {
            // [pid, name, [members]] -- keep member order as declared
            auto node = container["Container"];
            auto id = node[0].as<PID>();
            auto name = node[1].as<std::string>();
            std::vector<PID> members;
            for(const auto &member : node[2]) members.push_back(member.as<PID>());

            RegisterContainer(id, name, members);
        }
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

const std::vector<PID> &ParticleInfoEntry::Members() const noexcept {
    return members;
}

void ParticleInfoEntry::ClearMembers() noexcept {
    members.clear();
}

void ParticleInfoEntry::AddMember(PID pid) {
    if(std::find(members.begin(), members.end(), pid) != members.end()) return;

    // Mass-only container invariant: the first member sets the container's kinematics;
    // later members must agree in mass within tolerance.
    ParticleInfo part(pid);
    const double member_mass = part.Mass();
    if(members.empty()) {
        mass = member_mass;
        hmass = member_mass;
        massive = part.IsMassive();
    } else {
        const double group_mass = massive ? mass : 0.0;
        if(std::abs(member_mass - group_mass) > mass_tolerance) {
            auto msg = fmt::format(
                "ParticleInfo: adding particle of mass {} to container of mass {} exceeds "
                "tolerance {}",
                member_mass, group_mass, mass_tolerance);
            throw std::runtime_error(msg);
        }
    }
    members.push_back(pid);
}

ParticleInfo::ParticleDB ParticleInfo::particleDB;
std::map<std::string, achilles::PID> ParticleInfo::nameToPID;
std::map<PID, ParticleInfo::ContainerRule> ParticleInfo::rules;

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

ParticleInfo ParticleInfo::operator[](size_t i) const {
    if(!IsGroup()) return *this;
    if(i >= Size()) {
        auto msg =
            fmt::format("ParticleInfo: Asking for group member {} for group of size {}", i, Size());
        throw std::runtime_error(msg);
    }
    // Members are stored as particles; bar them when this container handle is anti
    ParticleInfo member(info->Members()[i]);
    return anti ? member.Anti() : member;
}

bool ParticleInfo::Includes(const ParticleInfo &other) const {
    // Compares by ID() (unsigned), so it does not distinguish a member from its anti-particle.
    return std::find(info->Members().begin(), info->Members().end(), other.ID()) !=
           info->Members().end();
}

std::vector<ParticleInfo> ParticleInfo::Decompose() const {
    if(!IsGroup()) return {*this};
    std::vector<ParticleInfo> parts;
    for(const auto &member : info->Members()) {
        ParticleInfo part(member);
        parts.push_back(anti ? part.Anti() : part);
    }
    return parts;
}

bool ParticleInfo::IsStable() const noexcept {
    if(info->stable == 0) return false;
    if(info->stable == 1) return true;
    if(info->stable == 2 && !IsAnti()) return true;
    if(info->stable == 3 && IsAnti()) return true;
    return false;
}

void ParticleInfo::RegisterContainer(PID id, std::string name, std::vector<PID> members) {
    auto entry = std::make_shared<ParticleInfoEntry>(ParticleInfoEntry());
    entry->id = id;
    entry->idname = name;
    entry->antiname = name;
    // AddMember adopts the mass and enforces the mass-only tolerance.
    for(const auto &member : members) entry->AddMember(member);

    particleDB[id] = entry;
    nameToPID[name] = id;
}

void ParticleInfo::RegisterContainerRule(PID id, std::string name, ContainerRule pred) {
    auto entry = std::make_shared<ParticleInfoEntry>(ParticleInfoEntry());
    entry->id = id;
    entry->idname = name;
    entry->antiname = name;

    particleDB[id] = entry;
    nameToPID[name] = id;
    rules[id] = pred;
}

void ParticleInfo::FinalizeContainers() {
    // Rebuild rule-based membership from the finalized table.
    for(const auto &[group_id, rule] : rules) {
        auto &group = particleDB[group_id];
        group->ClearMembers();
        for(const auto &[pid, entry] : particleDB) {
            if(entry->IsGroup() || rules.count(pid)) continue;
            if(rule(ParticleInfo(pid))) group->AddMember(pid);
        }
    }
}
