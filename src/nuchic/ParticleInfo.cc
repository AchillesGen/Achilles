#include "nuchic/ParticleInfo.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

#include <iostream>
#include <memory>

using nuchic::PID;
using nuchic::ParticleInfoEntry;
using nuchic::ParticleInfo;

// nuchic::Database::ParticleDB nuchic::Database::particleDB{};

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

namespace YAML {
template<>
struct convert<ParticleInfoEntry> {
    static Node encode(const ParticleInfoEntry &partInfo) {
        Node node;
        node.push_back(static_cast<int>(partInfo.id));
        node.push_back(partInfo.mass);
        node.push_back(partInfo.width);
        node.push_back(partInfo.icharge);
        node.push_back(partInfo.strong);
        node.push_back(partInfo.spin);
        node.push_back(partInfo.stable);
        node.push_back(partInfo.majorana);
        node.push_back(partInfo.massive);
        node.push_back(partInfo.hadron);
        node.push_back(partInfo.idname);
        node.push_back(partInfo.antiname);

        return node;
    }

    static bool decode(const Node &node, ParticleInfoEntry &partInfo) {
        if(!node.IsSequence() || node.size() != 12) {
            return false;
        } 

        partInfo.id = static_cast<nuchic::PID>(node[0].as<int>());
        partInfo.mass = node[1].as<double>();
        partInfo.width = node[2].as<double>();
        partInfo.icharge = node[3].as<int>();
        partInfo.strong = node[4].as<int>();
        partInfo.spin = node[5].as<int>();
        partInfo.stable = node[6].as<int>();
        partInfo.majorana = node[7].as<int>();
        partInfo.massive = node[8].as<bool>();
        partInfo.hadron = node[9].as<bool>();
        partInfo.idname = node[10].as<std::string>();
        partInfo.antiname = node[11].as<std::string>();

        return true;
    }
};
}

void nuchic::ParticleInfo::BuildDatabase(const std::string &datafile) {
    YAML::Node particleYAML = YAML::LoadFile(datafile);
    auto particles = particleYAML["Particles"];
    for(auto particle : particles) {
        auto entry = std::make_shared<ParticleInfoEntry>(particle["Particle"].as<ParticleInfoEntry>());
        particleDB.emplace(entry -> id, entry);
    }
    PrintDatabase();
}

void nuchic::ParticleInfo::PrintDatabase() {
    fmt::print("{:>7s} {:<20s} {:<20s} {:^10s}\t{:^10s}\n",
               "PID", "Name", "Anti-name", "Mass (MeV)", "Width (MeV)");
    for(const auto &part : particleDB) {
        std::cout <<  *(part.second) << "\n";
    }
}

ParticleInfoEntry::ParticleInfoEntry(const nuchic::PID &id_, const double &mass_,
                                     const double &width_, const int &icharge_,
                                     const int &strong_, const int &spin_,
                                     const int &stable_, const int &majorana_, 
                                     const bool &massive_, const bool &hadron_,
                                     std::string idname_, std::string antiname_) 
    : id(id_), mass(mass_), hmass(mass_), width(width_), icharge(icharge_), strong(strong_),
      spin(spin_), stable(stable_), majorana(majorana_), massive(massive_), hadron(hadron_), 
      idname(std::move(idname_)), antiname(std::move(antiname_)) {}

std::ostream& operator<<(std::ostream &os, const ParticleInfoEntry &entry) {
    os << entry.ToString();
    return os;
}

ParticleInfo::ParticleDB ParticleInfo::particleDB;

bool ParticleInfo::IsBaryon() const noexcept {
    if(IntID() % 10000 < 1000) return false;
    return true;
}

bool ParticleInfo::IsBHadron() const noexcept {
    if(IntID() < 100) return false;
    if(IntID()-100*IntID()/100 < 10) return false;
    if((IntID()-100*IntID()/100) / 10 == 5) return true;
    if((IntID()-1000*IntID()/1000) / 100 == 5) return true;
    if((IntID()-10000*IntID()/10000) / 1000 == 5) return true;
    return false;
}

bool ParticleInfo::IsCHadron() const noexcept {
    if(IntID() < 100) return false;
    if(IntID() - 100*IntID()/100 < 10) return false;
    if((IntID()-100*IntID()/100) / 10 == 4) return true;
    if((IntID()-1000*IntID()/1000) / 100 == 4) return true;
    if((IntID()-10000*IntID()/10000) / 1000 == 4) return true;
    return false;
}

double ParticleInfo::GenerateLifeTime() const {
    throw std::runtime_error("Not Implemented Yet");
    return 0.0;
}

bool ParticleInfo::IsStable() const noexcept {
    if(info -> stable == 0) return false;
    if(info -> stable == 1) return true;
    if(info -> stable == 2 && !IsAnti()) return true;
    if(info -> stable == 3 && IsAnti()) return true;
    return false;
}
