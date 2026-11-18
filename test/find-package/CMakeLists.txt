# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

# Downstream project used by CI to check that an installed Achilles works via find_package.
cmake_minimum_required(VERSION 3.17)
project(AchillesInstallTest LANGUAGES CXX)

find_package(Achilles REQUIRED)

add_executable(achilles-install-test main.cc)
target_link_libraries(achilles-install-test PRIVATE achilles::Achilles)
