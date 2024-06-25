#pragma once

#include <stdexcept>

namespace achilles {

class AchillesLoadError : public std::runtime_error {
  public:
    AchillesLoadError(const std::string &filename)
        : std::runtime_error("Achilles: Could not load " + filename) {}
};

} // namespace achilles
