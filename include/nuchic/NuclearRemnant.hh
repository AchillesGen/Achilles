#ifndef NUCLEAR_REMNANT_HH
#define NUCLEAR_REMNANT_HH

#include "spdlog/fmt/ostr.h"
#include "nuchic/Constants.hh"

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

        int PID() const { return std::stoi(fmt::format("100{:03}{:03}0", m_nZ, m_nA)); }
        double Mass() const { return static_cast<double>(m_nA)*Constant::mN; }

    private:
        size_t m_nA{}, m_nZ{};

};

}

#endif
