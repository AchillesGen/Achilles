# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(Manager)
add_subdirectory(Cuts)
# add_subdirectory(NucleonInteractions)

if(ACHILLES_ENABLE_SHERPA)
add_subdirectory(Sherpa)
endif()

if(ACHILLES_ENABLE_HEPMC3)
add_subdirectory(HepMC3)
add_subdirectory(NuHepMC)
endif()
