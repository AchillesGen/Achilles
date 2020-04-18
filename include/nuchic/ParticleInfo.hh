#ifndef PARTICLE_ID_HH
#define PARTICLE_ID_HH

// The classes in this file are inspired from the implementation found in Sherpa

#include <unordered_map>
#include <memory>
#include <spdlog/spdlog.h>
#include <string>
#include <functional>

namespace YAML {
    template<typename T>
    struct convert;
}

namespace nuchic {
    class ParticleInfoEntry;
    // class ParticleInfo;
}

std::ostream& operator<<(std::ostream&, const nuchic::ParticleInfoEntry&);
// std::ostream& operator<<(std::ostream&, const nuchic::ParticleInfo&);

namespace nuchic {
    class PID {
        public:
            constexpr PID(const long int &id_) : id(id_) {}
            // PID(const long int&, bool);
            constexpr long int AsInt() const { return id; }
            constexpr PID Anti() const { return PID(-id); }
            constexpr bool operator==(const PID &other) const { return id == other.id; }
            constexpr bool operator!=(const PID &other) const { return id != other.id; }
            constexpr bool operator<(const PID &other) const { return id < other.id; }
            constexpr bool operator>(const PID &other) const { return id > other.id; }
            constexpr operator long int() const { return id; }
            constexpr operator int() const { return static_cast<int>(id); }

            // Names for common particles
            // undefined
            static constexpr PID undefined() { return PID(0); }
            // Quarks
            static constexpr PID down() { return PID(1); } 
            static constexpr PID up() { return PID(2); } 
            static constexpr PID strange() { return PID(3); }
            static constexpr PID charm() { return PID(4); }
            static constexpr PID bottom() { return PID(5); }
            static constexpr PID top() { return PID(6); }
            // Leptons
            static constexpr PID electron() { return PID(11); }
            static constexpr PID nu_electron() { return PID(12); }
            static constexpr PID muon() { return PID(13); }
            static constexpr PID nu_muon() { return PID(14); }
            static constexpr PID tau() { return PID(15); }
            static constexpr PID nu_tau() { return PID(16); }
            // Gauge Bosons
            static constexpr PID gluon() { return PID(21); }
            static constexpr PID photon() { return PID(22); }
            static constexpr PID Zboson() { return PID(23); }
            static constexpr PID Wboson() { return PID(24); }
            static constexpr PID Higgs() { return PID(25); }
            // Mesons
            static constexpr PID pion0() { return PID(111); }
            static constexpr PID pionp() { return PID(211); }
            // Baryons
            static constexpr PID proton() { return PID(2212); }
            static constexpr PID neutron() { return PID(2112); }

        private:
            // Ensure id matches the numbering scheme defined at:
            // http://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
            // TODO: Implement this function, may not be explicitly needed though
            // bool valid(const long int&);
            long int id;
    };

    class ParticleInfo;

    class ParticleInfoEntry {
        public:
            ParticleInfoEntry()
                : id(PID::undefined()), mass(0.0), hmass(0.0), width(0.0), icharge(0),
                  strong(0), spin(0), stable(1), majorana(0), massive(false), hadron(false) {}
            ParticleInfoEntry(const PID&, const double&, const double&, const int&, const int&,
                              const int&, const int&, const int&, const bool&, const bool&,
                              const std::string&, const std::string&);

            ParticleInfoEntry(const ParticleInfoEntry &info) = default;
            ParticleInfoEntry(ParticleInfoEntry &&info) = default;

            ~ParticleInfoEntry() = default;

            bool operator==(const ParticleInfoEntry &other) {
                return id == other.id && mass == other.mass && hmass == other.hmass && width == other.width
                    && icharge == other.icharge && strong == other.strong && spin == other.spin && stable == other.stable
                    && majorana == other.majorana && massive == other.massive && hadron == other.hadron && idname == other.idname
                    && antiname == other.antiname;
            }

            friend class ParticleInfo;
            friend struct YAML::convert<ParticleInfoEntry>;
            friend std::ostream& ::operator<<(std::ostream&, const ParticleInfoEntry&);

        private:
            PID id;
            double mass, hmass, width;
            int icharge, strong;
            int spin, stable, majorana;
            bool massive;
            bool hadron;
            std::string idname, antiname;
    };

    struct PIDHasher {
        std::size_t operator()(const PID& key) const noexcept {
            return std::hash<long int>{}(key.AsInt()); 
        }
    };

    class ParticleInfo {
        private:
            using ParticleDB = std::unordered_map<PID, std::shared_ptr<ParticleInfoEntry>, PIDHasher>;
            static ParticleDB particleDB;
            static void BuildDatabase(const std::string&);

        public:
            ParticleInfo(const std::shared_ptr<ParticleInfoEntry> &info_, const bool &anti_=false)
                : info(info_), anti(false) {
                InitDatabase("data/Particles.yml");
                if(anti && info -> majorana == 0) anti = anti_;
            }

            ParticleInfo(const long int &id) : info(nullptr), anti(false) {
                InitDatabase("data/Particles.yml");
                auto it(particleDB.find(static_cast<PID>(std::abs(id))));
                if(it != particleDB.end()) 
                    info = it -> second;
                else
                    throw std::runtime_error(fmt::format("Invalid PID: id={}", id));
                if(id < 0 && info -> majorana == 0) anti = 1;
            }

            ParticleInfo(const PID &id, const bool &anti_) : info(nullptr), anti(false) {
                InitDatabase("data/Particles.yml");
                auto it(particleDB.find(id));
                info = it -> second;
                if(anti_ && info -> majorana == 0) anti = anti_;
            }

            ParticleInfo(const ParticleInfo&) = default;
            ParticleInfo(ParticleInfo&&) = default;

            ~ParticleInfo() = default;

            ParticleInfo operator=(const ParticleInfo &other) {
                if(*this != other) {
                    info = other.info;
                    anti = other.anti;
                }
                return *this;
            }

            // Property functions
            std::string IDName() const { return anti ? info -> antiname : info -> idname; }
            PID ID() const { return info -> id; }
            int IntID() const { return static_cast<int>(info ->id); }
            bool IsBaryon() const;
            bool IsHadron() const { return info -> hadron; }
            bool IsBHadron() const;
            bool IsCHadron() const;
            bool IsAnti() const { return anti; }
            bool IsFermion() const { return IntSpin()%2 == 1; }
            bool IsBoson() const { return IntSpin()%2 == 0; }
            bool IsScalar() const { return IntSpin() == 0; }
            bool IsVector() const { return IntSpin() == 2; }
            bool IsTensor() const { return IntSpin() == 4; }
            bool IsPhoton() const { return info -> id == PID::photon(); }
            bool IsLepton() const { return IntID() > 10 && IntID() < 19; }
            bool IsQuark() const { return IntID() < 10; }
            bool IsGluon() const { return info -> id == PID::gluon(); }

            int IntCharge() const { int charge(info -> icharge); return anti ? -charge : charge; }
            double Charge() const { return static_cast<double>(IntCharge()) / 3.0; }
            int IntSpin() const { return info -> spin; }
            double Spin() const { return static_cast<double>(info -> spin)/2.0; }
            bool SelfAnti() const { return info -> majorana != 0; } 
            bool Majorana() const { return info -> majorana == 1; }
            int Stable() const { return info -> stable; }
            bool IsStable() const;
            bool IsMassive() const { return info -> mass != 0 ? info -> massive : 0; } 
            double Mass() const { return info -> massive ? info -> mass : 0.0; }
            double Width() const { return info -> width; }

            double GenerateLifeTime() const;

            bool operator==(const ParticleInfo &other) const {
                return info == other.info && anti == other.anti;
            }
            bool operator!=(const ParticleInfo &other) const { return !(*this == other); }

            static ParticleDB Database() { return particleDB; }
            static void InitDatabase(const std::string &filename) {
                if(particleDB.size() == 0) {
                    particleDB[PID::undefined()] =  std::make_shared<ParticleInfoEntry>(ParticleInfoEntry());
                    BuildDatabase(filename);
                }
            }
            static void PrintDatabase();

        private:
            std::shared_ptr<ParticleInfoEntry> info;
            bool anti;
    };

    namespace Database {
        static constexpr auto Particle = ParticleInfo::Database;
    }

}

#endif
