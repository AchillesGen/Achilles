// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <stdexcept>

namespace achilles {

class InvalidChannel : public std::runtime_error {
  public:
    InvalidChannel(const std::string &msg) : std::runtime_error(msg) {}
};

} // namespace achilles
