#include "Achilles/CombinedCuts.hh"
#include "Achilles/Particle.hh"

bool achilles::CutCollection::EvaluateCuts(const std::vector<achilles::Particle> &parts) {
    ntot++;
    bool result = true;
    spdlog::trace("Evaluating Cuts");
    for(size_t i = 0; i < parts.size(); ++i) {
        if(!parts[i].IsFinal() && !parts[i].IsPropagating()) continue;
        spdlog::trace("Making cut for {}", parts[i].ID());
        // Single Particle Cuts
        for(const auto &cut : one_part_cuts)
            if(cut.Contains(parts[i].ID()))
                result &= cut.MakeCut(parts[i].Momentum());
        
        // Two Particle Cuts
        for(size_t j = i+1; j < parts.size(); ++j) {
            if(!parts[j].IsFinal() && !parts[i].IsPropagating()) continue;
            for(const auto &cut : two_part_cuts)
                if(cut.Contains(parts[i].ID(), parts[j].ID()))
                    result &= cut.MakeCut(parts[i].Momentum(), parts[j].Momentum());
        }
    }

    if(result) npass++;
    return result;
}

double achilles::CutCollection::CutEfficiency() const {
    return static_cast<double>(npass)/static_cast<double>(ntot);
}

bool achilles::CutCollection::AddCut(const std::set<PID> &pids, std::unique_ptr<OneParticleCut> cut) {
    bool combined = false;
    for(auto & combined_cut : one_part_cuts) {
        if(combined_cut.m_pids == pids) {
            combined_cut.m_cuts.push_back(std::move(cut));
            combined = true;
        }
    }
    if(!combined) {
        CombinedOneParticleCut combined_cut;
        combined_cut.m_pids = pids;
        combined_cut.m_cuts.push_back(std::move(cut));
        one_part_cuts.push_back(std::move(combined_cut));
    }

    return true;
}

bool achilles::CutCollection::AddCut(const std::set<PID> &pids, std::unique_ptr<TwoParticleCut> cut) {
    bool combined = false;
    for(auto & combined_cut : two_part_cuts) {
        if(combined_cut.m_pids == pids) {
            combined_cut.m_cuts.push_back(std::move(cut));
            combined = true;
        }
    }
    if(!combined) {
        CombinedTwoParticleCut combined_cut;
        combined_cut.m_pids = pids;
        combined_cut.m_cuts.push_back(std::move(cut));
        two_part_cuts.push_back(std::move(combined_cut));
    }
    return true;
}
