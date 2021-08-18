#ifndef CUTS_HH
#define CUTS_HH

#include <utility>

#include "nuchic/ParticleInfo.hh"
#include "fmt/ranges.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace nuchic {

class FourVector;

enum class CutType {
    Q2,
    Energy,
    Momentum,
    InvariantMass,
    TransverseMomentum,
    AngleTheta,
    AnglePhi,
    ETheta2
};

template<typename OStream>
OStream& operator<<(OStream &os, const CutType &cutType) {
    switch(cutType) {
        case CutType::Q2:
            os << "nuchic::CutType::Q2";
            break;
        case CutType::Energy:
            os << "nuchic::CutType::Energy";
            break;
        case CutType::Momentum:
            os << "nuchic::CutType::Momentum";
            break;
        case CutType::InvariantMass:
            os << "nuchic::CutType::InvariantMass";
            break;
        case CutType::TransverseMomentum:
            os << "nuchic::CutType::TransverseMomentum";
            break;
        case CutType::AngleTheta:
            os << "nuchic::CutType::AngleTheta";
            break;
        case CutType::AnglePhi:
            os << "nuchic::CutType::AnglePhi";
            break;
        case CutType::ETheta2:
            os << "nuchic::CutType::ETheta2";
            break;
    }
    return os;
}

template<class T>
class CutBase {
    public:
        using cut_range = std::pair<T, T>;
        using cut_ranges = std::vector<std::pair<T, T>>;
        CutBase() = default;

        CutBase(const YAML::Node &node) {
            if(node["min"] && node["max"] && node["range"]) {
                throw std::runtime_error("CutRange: Invalid syntax. Got min, max, and range");
            } else if(node["min"] || node["max"]) {
                T min = std::numeric_limits<T>::lowest();
                T max = std::numeric_limits<T>::max();
                if(node["min"]) {
                    if(!node["min"].IsScalar())
                        throw std::runtime_error("CutRange: Invalid min value");
                    min = node["min"].as<T>();
                }
                if(node["max"]) {
                    if(!node["max"].IsScalar())
                        throw std::runtime_error("CutRange: Invalid max value");
                    max = node["max"].as<T>();
                }
                m_range = {{min, max}}; 
            } else if(node["range"]) {
                m_range = node["range"].as<cut_ranges>();
            } else {
                throw std::runtime_error("CutRange: Invalid syntax. Missing cut values");
            }
            if(m_range.size() == 1) {
                spdlog::trace("Found cut range: [{}, {}]", m_range[0].first, m_range[0].second);
            } else {
                std::string ranges{};
                for(const auto &range : m_range) {
                    ranges += fmt::format("[{}, {}], ", range.first, range.second);
                }
                spdlog::trace("Found cut range: [{}]", ranges.substr(0, ranges.size()-2));
            }
        }

    protected:
        bool CheckCut(const T &val) const {
            bool result = false;
            // Run through possible cut ranges, keeping the event if a single cut range
            // is satisfied. Note that these cut ranges are all for a single variable,
            // e.g., theta in [0,1] OR [2,3] OR ...
            // Return true if any single cut range is satisfied.
            for(const auto &cut : m_range) {
                result |= (val > cut.first && val < cut.second);
            }

            return result;
        }

    private:
        cut_ranges m_range;
};

template<class Base>
class CutFactory {
    using Constructor = std::function<std::unique_ptr<Base>(const YAML::Node&)>;

    static std::map<std::string, Constructor>& Registry() {
        static std::map<std::string, Constructor> registry;
        return registry;
    }

    public:
        static std::unique_ptr<Base> InitializeCut(const std::string &name, const YAML::Node &node) {
            auto constructor = Registry().at(name);
            return constructor(node);
        }

        template<class Derived>
        static void Register(std::string name) {
            if(Registry().find(name) != Registry().end())
                spdlog::error("{} is already registered!", name);
            spdlog::trace("Registering {} Cut", name);
            Registry()[name] = Derived::Construct;
        }

        static bool IsRegistered(std::string name) {
            return Registry().find(name) != Registry().end();
        }

        static void DisplayCuts() {
            fmt::print("Registered cuts:\n");
            for(auto &registered : Registry()) 
                fmt::print("  - {}\n", registered.first);
        }
};

template<class Base, class Derived>
class RegistrableCut {
    protected:
        RegistrableCut() = default;
        virtual ~RegistrableCut()  {
            if(!m_registered)
                spdlog::error("Error registering cut");
        }

        static bool Register() {
            CutFactory<Base>::template Register<Derived>(Derived::Name());
            return true;
        }

    private:
        static const bool m_registered;
};
template<class Base, class Derived>
const bool RegistrableCut<Base, Derived>::m_registered = RegistrableCut<Base, Derived>::Register();

class Cut {
    using cut_range = std::vector<std::pair<double, double>>;

    public:
        Cut() = default;

        // Calculate all cuts for FourVector defined by set of PIDs
        bool operator()(const FourVector&) const;

        // Add new cuts
        void Add(const Cut &other) {
            for(const auto &cut : other.m_cuts) {
                Add(cut.first, cut.second); 
            }
        }
        void Add(const CutType&, double);
        void Add(const CutType&, cut_range);
        void Add(const std::function<double(const FourVector&)>&, double);
        void Add(const std::function<double(const FourVector&)>&,
                 cut_range);

    private:
        bool MakeCut(double, const cut_range&) const;
        bool Q2Cut(const FourVector&, const cut_range&) const;
        bool EnergyCut(const FourVector&, const cut_range&) const;
        bool MomentumCut(const FourVector&, const cut_range&) const;
        bool InvariantMassCut(const FourVector&, const cut_range&) const;
        bool TransverseMomentumCut(const FourVector&, const cut_range&) const;
        bool AngleThetaCut(const FourVector&, const cut_range&) const;
        bool AnglePhiCut(const FourVector&, const cut_range&) const;
        bool ETheta2Cut(const FourVector&, const cut_range&) const;

        std::map<CutType, cut_range> m_cuts;
};

using Cuts = std::map<PID, Cut>;

}

namespace YAML {
template<>
struct convert<nuchic::CutType> {
    static bool decode(const Node &node, nuchic::CutType &cut) {
        auto name = node.as<std::string>();
        if(name == "Q2")
            cut = nuchic::CutType::Q2;
        else if(name == "Energy")
            cut = nuchic::CutType::Energy;
        else if(name == "Momentum")
            cut = nuchic::CutType::Momentum;
        else if(name == "Mass")
            cut = nuchic::CutType::InvariantMass;
        else if(name == "TransverseMomentum")
            cut = nuchic::CutType::TransverseMomentum;
        else if(name == "Theta")
            cut = nuchic::CutType::AngleTheta;
        else if(name == "Phi")
            cut = nuchic::CutType::AnglePhi;
        else if(name == "ETheta2")
            cut = nuchic::CutType::ETheta2;
        else
            return false;

        return true;
    }
};

template<typename T>
struct convert<std::pair<T, T>> {
    static bool decode(const Node &node, std::pair<T, T> &range) {
        if(node[0].IsScalar() && node[1].IsScalar() && node.size() == 2) {
            range = {node[0].as<T>(), node[1].as<T>()};
            return true;
        }
        return false;
    }
};

template<typename T>
struct convert<std::vector<std::pair<T, T>>> {
    static bool decode(const Node &node, std::vector<std::pair<T, T>> &cutRange) {
        if(!node.IsSequence())
            return false;
        if(node[0].IsScalar() && node[1].IsScalar() && node.size() == 2) {
            cutRange = {{node[0].as<double>(), node[1].as<T>()}};
        } else {
            for(const auto &subNode : node) {
                cutRange.push_back(subNode.as<std::pair<T, T>>());
            }
        }
        return true;
    }
};

template<>
struct convert<nuchic::Cut> {
    static bool decode(const Node &node, nuchic::Cut &cut) {
        auto type = node["type"].as<nuchic::CutType>();
        if(node["min"] && node["max"] && node["range"]) {
            fmt::format("Allowed options are a range, a min value, a max value, or both min and max values, not all three\n");
            return false;
        } else if(node["min"] && node["max"]) {
            if(!(node["min"].IsScalar() && node["max"].IsScalar()))
                return false;
            const auto min = node["min"].as<double>();
            const auto max = node["max"].as<double>();
            cut.Add(type, {{min, max}});
        } else if(node["min"]) {
            if(!node["min"].IsScalar())
                return false;
            const auto min = node["min"].as<double>();
            cut.Add(type, min);
        } else if(node["max"]) {
            if(!node["max"].IsScalar())
                return false;
            const auto max = node["max"].as<double>();
            cut.Add(type, {{-std::numeric_limits<double>::infinity(), max}});
        } else if(node["range"]) {
            const auto range = node["range"].as<std::vector<std::pair<double, double>>>();
            cut.Add(type, range);
        } else {
            fmt::format("Missing either a range, a min value, a max value, or both min and max values\n");
            return false;
        }

        return true;
    }
};

template<>
struct convert<nuchic::Cuts> {
    static bool decode(const Node &node, nuchic::Cuts &cuts) {
        for(const auto &subNode : node) {
            nuchic::PID pid;
            try{
                pid = nuchic::PID(subNode.first.as<int>());
            } catch(const std::exception &e) {
                auto name = subNode.first.as<std::string>();
                pid = nuchic::ParticleInfo::NameToPID().at(name); 
            }
            cuts[pid] = nuchic::Cut();
            for(const auto &cutNode : subNode.second) {
                cuts[pid].Add(cutNode["cut"].as<nuchic::Cut>());
            }
        }
        return true;
    }
};

}

#endif
