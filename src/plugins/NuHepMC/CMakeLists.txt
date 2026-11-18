# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(AchillesNuHepMC SHARED NuHepMCWriter.cc)
target_link_libraries(AchillesNuHepMC PRIVATE achilles_options
                                      PUBLIC NuHepMC::CPPUtils fmt::fmt spdlog::spdlog yaml-cpp::yaml-cpp AchillesUtilities AchillesPhysics AchillesNuclearModels)
