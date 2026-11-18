# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(plugin_manager PluginManager.cc)
target_link_libraries(plugin_manager PRIVATE achilles_options achilles_warnings
                                     PUBLIC event_gen)
