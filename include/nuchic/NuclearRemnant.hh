#ifndef NUCLEAR_REMNANT_HH
#define NUCLEAR_REMNANT_HH

#include "spdlog/fmt/ostr.h"

namespace nuchic {

class NuclearRemnant {
    public:
        NuclearRemnant() = default;
        NuclearRemnant(size_t nA, size_t nZ) : m_nA{std::move(nA)}, m_nZ{std::move(nZ)} {}

        template<typename OStream>
        friend OStream& operator<<(OStream &os, const NuclearRemnant &nuclearRemnant) { 
            os << "NuclearRemnant(" << nuclearRemnant.m_nA << ", " << nuclearRemnant.m_nZ << ")";
            return os;
        }

    private:
        size_t m_nA{}, m_nZ{};

};

}

#endif
