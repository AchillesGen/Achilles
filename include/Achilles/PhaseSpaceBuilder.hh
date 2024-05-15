#ifndef PHASE_SPACE_BUILDER
#define PHASE_SPACE_BUILDER

#include "Achilles/Achilles.hh"
#include "Achilles/PhaseSpaceMapper.hh"
#include "Achilles/ProcessInfo.hh"

#include <optional>

namespace PHASIC {
class Channels;
}

namespace achilles {

class Beam;

class PSBuilder {
  public:
    PSBuilder(const ProcessInfo &);
    MOCK ~PSBuilder() = default;
    MOCK PSBuilder &Beam(std::shared_ptr<Beam>, size_t = 0);
    MOCK PSBuilder &Hadron(const std::string &, size_t = 1);
    MOCK PSBuilder &FinalState(const std::string &, std::optional<double> = std::nullopt);
#ifdef ACHILLES_SHERPA_INTERFACE
    MOCK PSBuilder &SherpaFinalState(const std::string &);
    MOCK PSBuilder &GenFinalState(std::unique_ptr<PHASIC::Channels>);
#endif // ACHILLES_SHERPA_INERFACE

    MOCK std::unique_ptr<PSMapper> build() { return std::move(phase_space); }

  private:
    ProcessInfo m_info;
    [[maybe_unused]] size_t m_nlep, m_nhad, m_nspec;
    std::unique_ptr<PSMapper> phase_space = nullptr;
};

} // namespace achilles

#endif
