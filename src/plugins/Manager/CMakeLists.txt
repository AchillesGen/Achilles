# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(AchillesManager PluginManager.cc)
target_link_libraries(AchillesManager PRIVATE achilles_options achilles_warnings
                                      PUBLIC AchillesEventGen)
