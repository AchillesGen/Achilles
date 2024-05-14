#pragma once

#include "ATOOLS/Math/Vector.H"
#include "spdlog/fmt/ostr.h"
#include "fmt/core.h"

template <typename OStream> OStream &operator<<(OStream &os, const ATOOLS::Vec4D &vec4) {
    os << "Vec4D(" << vec4[0] << ", " << vec4[1] << ", " << vec4[2] << ", " << vec4[3] << ")";
    return os;
}

template<> struct fmt::formatter<ATOOLS::Vec4D> {
    char presentation = 'e';
    constexpr auto parse(format_parse_context& ctx) -> format_parse_context::iterator {
        // Parse the presentation format and store it in the formatter:
        auto it = ctx.begin(), end = ctx.end();
        if(it != end && (*it == 'f' || *it == 'e')) presentation = *it++;

        // Check if reached the end of the range:
        if(it != end && *it != '}')
            throw format_error("Invalid format");

        // Return an iterator past the end of the parsed range:
        return it;
    }

    auto format(const ATOOLS::Vec4D& p, format_context& ctx) const -> format_context::iterator {
        // ctx.out() is an output iterator to write to
        return format_to(
                ctx.out(),
                presentation == 'f' ? "ATOOLS::Vec4D({:.8f}, {:.8f}, {:.8f}, {:.8f})"
                                    : "ATOOLS::Vec4D({:.8e}, {:.8e}, {:.8e}, {:.8e})",
                p[0], p[1], p[2], p[3]);
    }
};
