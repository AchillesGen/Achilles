#ifndef COMBINED_CUTS_HH
#define COMBINED_CUTS_HH

#include "Achilles/OneParticleCuts.hh"
#include "Achilles/TwoParticleCuts.hh"

namespace achilles {

class Particle;
class CutCollection;

class CombinedOneParticleCut {
  public:
    bool Contains(PID pid) const { return m_pids.find(pid) != m_pids.end(); }
    bool MakeCut(const FourVector &mom) const {
        bool result = true;
        for(const auto &cut : m_cuts) result &= cut->MakeCut(mom);
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
        return (m_pids.find(pid1) != m_pids.end()) && (m_pids.find(pid2) != m_pids.end());
    }
    bool MakeCut(const FourVector &mom1, const FourVector &mom2) const {
        bool result = true;
        for(const auto &cut : m_cuts) { result &= cut->MakeCut(mom1, mom2); }
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
    bool EvaluateCuts(const std::vector<Particle> &);
    double CutEfficiency() const;
    bool AddCut(const std::set<PID> &, std::unique_ptr<OneParticleCut>);
    bool AddCut(const std::set<PID> &, std::unique_ptr<TwoParticleCut>);

    friend YAML::convert<CutCollection>;

  private:
    size_t npass{}, ntot{};
    std::vector<CombinedOneParticleCut> one_part_cuts;
    std::vector<CombinedTwoParticleCut> two_part_cuts;
};

} // namespace achilles

namespace YAML {

template <> struct convert<achilles::CutCollection> {
    static bool decode(const Node &node, achilles::CutCollection &cuts) {
        spdlog::trace("Loading cuts");
        for(const auto &subnode : node) {
            auto cut_type = subnode["Type"].as<std::string>();
            auto pid = subnode["PIDs"].as<std::set<achilles::PID>>();
            spdlog::trace("Found cut: {}", cut_type);
            spdlog::trace("Found PIDs: [{}]", std::vector<int>(pid.begin(), pid.end()));
            if(achilles::CutFactory<achilles::OneParticleCut>::IsRegistered(cut_type)) {
                auto cut = achilles::CutFactory<achilles::OneParticleCut>::InitializeCut(cut_type,
                                                                                         subnode);
                cuts.AddCut(pid, std::move(cut));
            } else if(achilles::CutFactory<achilles::TwoParticleCut>::IsRegistered(cut_type)) {
                auto cut = achilles::CutFactory<achilles::TwoParticleCut>::InitializeCut(cut_type,
                                                                                         subnode);
                cuts.AddCut(pid, std::move(cut));
            } else {
                return false;
            }
        }
        return true;
    }
};

} // namespace YAML

#endif
