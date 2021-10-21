#ifndef PHASE_SPACE_BUILDER
#define PHASE_SPACE_BUILDER

#include "nuchic/PhaseSpaceMapper.hh"

namespace nuchic {

class Beam;

class PSBuilder {
    public:
        PSBuilder(size_t nlep=2, size_t nhad=2) 
            : m_nlep{nlep}, m_nhad{nhad} { phase_space = std::make_unique<PSMapper>(nlep, nhad); }
        PSBuilder& Beam(std::shared_ptr<Beam>, size_t=1);
        PSBuilder& Hadron(const std::string&, size_t=0);
        PSBuilder& FinalState(const std::string&, const std::vector<double>&);
        PSBuilder& SherpaFinalState(const std::string&, const std::vector<double>&);

        std::unique_ptr<PSMapper> build() { return std::move(phase_space); }

    private:
        size_t m_nlep, m_nhad;
        std::unique_ptr<PSMapper> phase_space = nullptr;
};

}

#endif
