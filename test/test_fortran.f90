program fortran_test
  use unit_test
  use test_spectral_interface
  use test_map_interface

  implicit none

  type(test_suite_type) :: test_suite

  ! example with specific suite
  call test_suite_init('Achilles Fortran Test Suite', test_suite)
  ! Test all test cases
  call test_spectral(test_suite)
  call test_map(test_suite)
  ! report the complete suite
  call test_suite_report(test_suite)
  ! finalize
  call test_suite_final(test_suite)

end program fortran_test
