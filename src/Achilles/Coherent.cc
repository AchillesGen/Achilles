#ifdef ACHILLES_EVENT_DETAILS
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE
#endif

#include "Achilles/Event.hh"
#include "Achilles/HardScattering.hh"
#include "Achilles/HardScatteringFactory.hh"
#include "Achilles/Particle.hh"

REGISTER_HARDSCATTERING(achilles::FCoherent);

achilles::HCurrents achilles::FCoherent::HadronicCurrents(Event &event, const FFInfoMap &,
                                                          const FFInfoMap &,
                                                          const FFInfoMap &coherentFF) const {
    auto pIn = event.Momentum().front();
    auto pOut = event.Momentum().back();
    auto qVec = event.Momentum()[1];
    for(size_t i = 2; i < event.Momentum().size() - 1; ++i) { qVec -= event.Momentum()[i]; }

    auto ffVals = EvalFormFactor(-qVec.M2() / 1.0_GeV / 1.0_GeV);

    // Calculate coherent contributions
    HCurrents results;
    for(const auto &formFactor : coherentFF) {
        auto ffVal = CouplingsFF(ffVals, formFactor.second);
        SPDLOG_TRACE("fcoh = {}", ffVal[3]);

        Current current;
        std::vector<std::complex<double>> subcur(4);
        for(size_t i = 0; i < 4; ++i) { subcur[i] = (pIn[i] + pOut[i]) * ffVal[3]; };
        current.push_back(subcur);
        results[2][formFactor.first] = current;
    }

    return results;
}
