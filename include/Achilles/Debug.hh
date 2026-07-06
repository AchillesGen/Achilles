// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

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
