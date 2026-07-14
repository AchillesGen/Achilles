# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(nuhepmc SHARED NuHepMCWriter.cc)
target_link_libraries(nuhepmc PRIVATE achilles_options achilles_warnings
                               PUBLIC NuHepMC::CPPUtils fmt::fmt spdlog::spdlog yaml::cpp utilities physics nuclear_models)
