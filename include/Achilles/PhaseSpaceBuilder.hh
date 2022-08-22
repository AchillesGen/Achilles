#ifndef PHASE_SPACE_BUILDER
#define PHASE_SPACE_BUILDER

#include "Achilles/Achilles.hh"
#include "Achilles/PhaseSpaceMapper.hh"

namespace PHASIC {
    class Channels;
}

namespace achilles {

class Beam;

class PSBuilder {
    public:
#ifdef ENABLE_BSM
        PSBuilder(size_t nlep=2, size_t nhad=2) 
            : m_nlep{nlep}, m_nhad{nhad} { phase_space = std::make_unique<PSMapper>(nlep, nhad); }
        MOCK PSBuilder& SherpaFinalState(const std::string&, const std::vector<double>&);
        MOCK PSBuilder& GenFinalState(std::unique_ptr<PHASIC::Channels>);
#else
        PSBuilder(size_t nlep=2, size_t nhad=2) { phase_space = std::make_unique<PSMapper>(nlep, nhad); }
#endif // ENABLE_BSM
        MOCK ~PSBuilder() = default;
        MOCK PSBuilder& Beam(std::shared_ptr<Beam>, const std::vector<double>&, size_t=1);
        MOCK PSBuilder& Hadron(const std::string&, const std::vector<double>&, size_t=0);
        MOCK PSBuilder& FinalState(const std::string&, const std::vector<double>&);

        MOCK std::unique_ptr<PSMapper> build() { return std::move(phase_space); }

    private:
#ifdef ENABLE_BSM
        size_t m_nlep, m_nhad;
#endif
        std::unique_ptr<PSMapper> phase_space = nullptr;
};

}

#endif
