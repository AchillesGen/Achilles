! SPDX-FileCopyrightText: 2018-2026 Achilles Developers
! SPDX-License-Identifier: GPL-3.0-or-later

program fortran_test
  use unit_test
  use test_spectral_interface
  use test_map_interface
  use test_fortran_inteference
  use test_fortran_utils

  implicit none

  type(test_suite_type) :: test_suite
  logical, allocatable :: results(:)

  ! example with specific suite
  call test_suite_init('Achilles Fortran Test Suite', test_suite)
  ! Test all test cases
  call test_spectral(test_suite)
  call test_map(test_suite)
  call test_interference(test_suite)
  call test_string_conversion(test_suite)
  ! report the complete suite
  call test_suite_report(test_suite)
  ! capture per-assertion pass/fail results before the suite is torn down
  results = test_suite_get_assert_results(test_suite)
  ! finalize
  call test_suite_final(test_suite)

  ! The FUT framework always exits 0, so signal failures with a non-zero exit
  ! code. This lets CTest (and CI) detect a failing Fortran assertion.
  if (any(.not. results)) error stop 1

end program fortran_test
