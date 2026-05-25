// SPDX-FileCopyrightText: 2018-2024 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <stdexcept>
#include <string>

namespace achilles {

class AchillesLoadError : public std::runtime_error {
  public:
    AchillesLoadError(const std::string &filename)
        : std::runtime_error("Achilles: Could not load " + filename) {}
};

class AchillesCascadeError : public std::runtime_error {
  public:
    AchillesCascadeError(const std::string &msg) : std::runtime_error(msg) {}
};

} // namespace achilles
