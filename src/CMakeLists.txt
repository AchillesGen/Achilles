# SPDX-FileCopyrightText: 2018-2024 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

if(ACHILLES_ENABLE_GZIP)
add_subdirectory(gzstream)
endif()
add_subdirectory(plugins)
add_subdirectory(Achilles)
