#include "nuchic/CombinedCuts.hh"
#include "nuchic/Particle.hh"

bool nuchic::CutCollection::EvaluateCuts(const std::vector<nuchic::Particle> &parts) {
    ntot++;
    bool result = true;
    spdlog::trace("Evaluating Cuts");
    for(size_t i = 0; i < parts.size(); ++i) {
        spdlog::trace("Making cut for {}, status = {}", parts[i].ID(), parts[i].Status());
        if(!parts[i].IsFinal()) continue;
        spdlog::trace("Making cut for {}", parts[i].ID());
        // Single Particle Cuts
        for(const auto &cut : one_part_cuts)
            if(cut.Contains(parts[i].ID()))
                result &= cut.MakeCut(parts[i].Momentum());
        
        // Two Particle Cuts
        for(size_t j = i+1; j < parts.size(); ++j) {
            if(!parts[j].IsFinal()) continue;
            for(const auto &cut : two_part_cuts)
                if(cut.Contains(parts[i].ID(), parts[j].ID()))
                    result &= cut.MakeCut(parts[i].Momentum(), parts[j].Momentum());
        }
    }

    if(result) npass++;
    return result;
}

double nuchic::CutCollection::CutEfficiency() const {
    return static_cast<double>(npass)/static_cast<double>(ntot);
}

bool nuchic::CutCollection::AddCut(const std::set<PID> &pids, std::unique_ptr<OneParticleCut> cut) {
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

bool nuchic::CutCollection::AddCut(const std::set<PID> &pids, std::unique_ptr<TwoParticleCut> cut) {
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
