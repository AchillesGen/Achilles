#ifndef PRINT_VEC_HH
#define PRINT_VEC_HH

#include "ATOOLS/Math/Vector.H"
#include "spdlog/fmt/ostr.h"

template <typename OStream> OStream &operator<<(OStream &os, const ATOOLS::Vec4D &vec4) {
    os << "Vec4D(" << vec4[0] << ", " << vec4[1] << ", " << vec4[2] << ", " << vec4[3] << ")";
    return os;
}

#endif
