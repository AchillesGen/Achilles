// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SPLASH_HH
#define SPLASH_HH

#include <string>
#include <vector>

namespace achilles {

std::vector<std::string> SplashLines(const std::string &version);
void PrintSplash(const std::string &version);

} // namespace achilles

#endif
