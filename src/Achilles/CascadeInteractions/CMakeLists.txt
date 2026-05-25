# SPDX-FileCopyrightText: 2018-2024 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(AchillesCascadeInteractions SHARED
    GeantInteractions.cc
    NasaInteractions.cc
    ConstantInteractions.cc
    MesonBaryonInteractions.cc
    OsetMesonBaryonInteractions.cc
    PionAbsorption.cc
    PionInteractions.cc
    DeltaInteractions.cc
    NucleonNucleon.cc
)
target_link_libraries(AchillesCascadeInteractions PRIVATE project_options project_warnings
                                                  PUBLIC physics AchillesIntegrator)
list(APPEND achilles_targets AchillesCascadeInteractions)
set(achilles_targets ${achilles_targets} PARENT_SCOPE)
