#pragma once

#include "Achilles/FourVector.hh"

#include <string>
#include <vector>

namespace achilles {

struct DebugEvents {
    DebugEvents(const std::string &filename, size_t multiplicity);
    std::vector<std::vector<FourVector>> events;
};

} // namespace achilles
