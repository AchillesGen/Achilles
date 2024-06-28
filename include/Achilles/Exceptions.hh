#pragma once

#include <stdexcept>

namespace achilles {

class InvalidChannel : public std::runtime_error {
  public:
    InvalidChannel(const std::string &msg) : std::runtime_error(msg) {}
};

} // namespace achilles
