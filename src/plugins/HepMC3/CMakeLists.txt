# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(AchillesHepMC3 SHARED HepMC3EventWriter.cc)
target_link_libraries(AchillesHepMC3 PRIVATE achilles_options
                             PUBLIC HepMC3::HepMC3 fmt::fmt spdlog::spdlog yaml-cpp::yaml-cpp utilities physics nuclear_models)
