#include "Achilles/CombinedCuts.hh"
#include "Achilles/Particle.hh"

bool achilles::CutCollection::EvaluateCuts(const refParticles& parts) {
    ntot++;
    bool result = true;
    spdlog::trace("Evaluating Cuts");
    for(size_t i = 0; i < parts.size(); ++i) {
		Particle& part=parts[i].get();
        if(!part.IsFinal() && !part.IsPropagating()) continue;
		PID id=part.ID();
        spdlog::trace("Making cut for {}", id);
        // Single Particle Cuts
        for(const auto &cut : one_part_cuts)
            if(cut.Contains(id)) result &= cut.MakeCut(part.Momentum());

        // Two Particle Cuts
        for(size_t j = i + 1; j < parts.size(); ++j) {
			Particle& part_j=parts[j].get();
            if(!part_j.IsFinal() && !part.IsPropagating()) continue;
            for(const auto &cut : two_part_cuts)
                if(cut.Contains(part.ID(), part_j.ID()))
                    result &= cut.MakeCut(part.Momentum(), part_j.Momentum());
        }
    }

    if(result) npass++;
    return result;
}

double achilles::CutCollection::CutEfficiency() const {
    return static_cast<double>(npass) / static_cast<double>(ntot);
}

bool achilles::CutCollection::AddCut(const std::set<PID> &pids,
                                     std::unique_ptr<OneParticleCut> cut) {
    bool combined = false;
    for(auto &combined_cut : one_part_cuts) {
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

bool achilles::CutCollection::AddCut(const std::set<PID> &pids,
                                     std::unique_ptr<TwoParticleCut> cut) {
    bool combined = false;
    for(auto &combined_cut : two_part_cuts) {
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
