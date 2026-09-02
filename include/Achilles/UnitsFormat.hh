// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later
//
// fmt/spdlog formatting for achilles::units::Quantity (single-axis design).
//
//   Energy beam = 2.2_GeV;   Length r = 1.2_fm;
//   fmt::print("{}\n",     beam);          // "2200 MeV"
//   fmt::print("{}\n",     r);             // "0.00608 MeV^-1"  (canonical!)
//   fmt::print("{:.3f}\n", in(r,   fm));   // "1.200 fm"
//   fmt::print("{:.3f}\n", in(sig, nb));   // "11.681 nb"
//
// Default prints the canonical storage unit (always a power of MeV), built at
// compile time from the mass dimension D. A length prints as MeV^-1 -- that is
// how it is stored -- so use in(q, fm) when you want fm.
//
// spdlog is built with SPDLOG_FMT_EXTERNAL (see external/CMakeLists.txt), so
// these formatters are the ones spdlog sees.

#ifndef UNITS_FORMAT_HH
#define UNITS_FORMAT_HH

#include "Achilles/Units.hh"

#include <fmt/format.h>
#include <string_view>

// Default formatter: <number> MeV^D (D omitted when 1, unit omitted when 0).
template <int D> struct fmt::formatter<achilles::units::Quantity<D>> {
    fmt::formatter<double> num;

    constexpr auto parse(fmt::format_parse_context &ctx) { return num.parse(ctx); }

    template <typename Ctx> auto format(const achilles::units::Quantity<D> &q, Ctx &ctx) const {
        auto out = num.format(q.native(), ctx);
        if constexpr(D != 0) {
            out = fmt::format_to(out, " MeV");
            if constexpr(D != 1) out = fmt::format_to(out, "^{}", D);
        }
        return out;
    }
};

// Choose a display unit: in(q, unit) -> value already in that unit + its label.
// Dimension-checked by Quantity::in(): in(energy, fm) is a compile error.
namespace achilles::units {

template <int D> struct Scaled {
    double value;
    std::string_view label;
};

template <int D> constexpr Scaled<D> in(Quantity<D> q, Unit<D> u) noexcept {
    return Scaled<D>{q.in(u), u.name()};
}

} // namespace achilles::units

template <int D> struct fmt::formatter<achilles::units::Scaled<D>> {
    fmt::formatter<double> num;

    constexpr auto parse(fmt::format_parse_context &ctx) { return num.parse(ctx); }

    template <typename Ctx> auto format(const achilles::units::Scaled<D> &s, Ctx &ctx) const {
        auto out = num.format(s.value, ctx);
        if(!s.label.empty()) {
            *out++ = ' ';
            out = fmt::format_to(out, "{}", s.label);
        }
        return out;
    }
};

#endif
