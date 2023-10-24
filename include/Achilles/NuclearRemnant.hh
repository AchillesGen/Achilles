#ifndef NUCLEAR_REMNANT_HH
#define NUCLEAR_REMNANT_HH

#include "Achilles/Constants.hh"
#include "spdlog/fmt/ostr.h"

namespace achilles {

class NuclearRemnant {
  public:
    NuclearRemnant() = default;
    NuclearRemnant(size_t nA, size_t nZ) : m_nA{std::move(nA)}, m_nZ{std::move(nZ)} {}

    template <typename OStream>
    friend OStream &operator<<(OStream &os, const NuclearRemnant &nuclearRemnant) {
        os << "NuclearRemnant(" << nuclearRemnant.m_nA << ", " << nuclearRemnant.m_nZ << ")";
        return os;
    }

    int PID() const { return std::stoi(fmt::format("100{:03}{:03}0", m_nZ, m_nA)); }
    double Mass() const { return static_cast<double>(m_nA) * Constant::mN; }
    bool operator==(const NuclearRemnant &other) const {
        return m_nA == other.m_nA && m_nZ == other.m_nZ;
    }

  private:
    size_t m_nA{}, m_nZ{};
};

} // namespace achilles

namespace fmt {

template <> struct formatter<achilles::NuclearRemnant> {
    template <typename ParseContext> constexpr auto parse(ParseContext &ctx) { return ctx.begin(); }

    template <typename FormatContext> auto format(const achilles::NuclearRemnant &nucrem, FormatContext &ctx) const {
        std::stringstream ss;
        ss << nucrem;

        return format_to(ctx.out(), ss.str());
    }
};

} // namespace fmt

#endif
