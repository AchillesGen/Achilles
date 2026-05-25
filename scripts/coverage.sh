#!/bin/bash
# SPDX-FileCopyrightText: 2018-2024 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later


./test/achilles-testsuite
./test/achilles-fortran-testsuite

lcov --directory . --capture --output-file coverage.info --ignore-errors gcov,gcov --ignore-errors mismatch
lcov --remove coverage.info '/usr/*' 'build/*' 'test/*' --output-file coverage.info
genhtml coverage.info --output-directory coverage
