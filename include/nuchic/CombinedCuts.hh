#ifndef COMBINED_CUTS_HH
#define COMBINED_CUTS_HH

#include "nuchic/OneParticleCuts.hh"
#include "nuchic/TwoParticleCuts.hh"

namespace nuchic {

class Particle;
class CutCollection;

class CombinedOneParticleCut {
    public:
        bool Contains(PID pid) const { return m_pids.find(pid) != m_pids.end(); }
        bool MakeCut(const FourVector &mom) const {
            bool result = true;
            for(const auto &cut : m_cuts) result &= cut -> MakeCut(mom); 
            return result;
        }

        friend YAML::convert<CutCollection>;
        friend CutCollection;
    private:
        std::set<PID> m_pids;
        std::vector<std::unique_ptr<OneParticleCut>> m_cuts;
};

class CombinedTwoParticleCut {
    public:
        bool Contains(PID pid1, PID pid2) const {
            return (m_pids.find(pid1) != m_pids.end())
                && (m_pids.find(pid2) != m_pids.end());
        }
        bool MakeCut(const FourVector &mom1, const FourVector &mom2) const {
            bool result = true;
            for(const auto &cut : m_cuts) {
                result &= cut -> MakeCut(mom1, mom2);
            }
            return result;
        }

        friend YAML::convert<CutCollection>;
        friend CutCollection;
    private:
        std::set<PID> m_pids;
        std::vector<std::unique_ptr<TwoParticleCut>> m_cuts;
};

class CutCollection {
    public:
        bool EvaluateCuts(const std::vector<Particle>&);
        double CutEfficiency() const;
        bool AddCut(const std::set<PID>&, std::unique_ptr<OneParticleCut>);
        bool AddCut(const std::set<PID>&, std::unique_ptr<TwoParticleCut>);

        friend YAML::convert<CutCollection>;
    private:
        size_t npass{}, ntot{};
        std::vector<CombinedOneParticleCut> one_part_cuts;
        std::vector<CombinedTwoParticleCut> two_part_cuts;
};

}

namespace YAML {

template<>
struct convert<nuchic::CutCollection> {
    static bool decode(const Node &node, nuchic::CutCollection &cuts) {
        spdlog::trace("Loading cuts");
        for(const auto &subnode : node) {
            auto cut_type = subnode["Type"].as<std::string>();
            auto pid = subnode["PIDs"].as<std::set<nuchic::PID>>();
            spdlog::trace("Found cut: {}", cut_type);
            spdlog::trace("Found PIDs: [{}]", std::vector<int>(pid.begin(), pid.end()));
            if(nuchic::CutFactory<nuchic::OneParticleCut>::IsRegistered(cut_type)) {
                auto cut = nuchic::CutFactory<nuchic::OneParticleCut>::InitializeCut(cut_type, subnode);
                cuts.AddCut(pid, std::move(cut));
            } else if(nuchic::CutFactory<nuchic::TwoParticleCut>::IsRegistered(cut_type)) {
                auto cut = nuchic::CutFactory<nuchic::TwoParticleCut>::InitializeCut(cut_type, subnode);
                cuts.AddCut(pid, std::move(cut));
            } else {
                return false;
            }
        }
        return true;
    }
};

}

#endif
