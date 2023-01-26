program fortran_test
  use unit_test
  use test_spectral_interface

  implicit none

  type(test_suite_type) :: test_suite

  ! example with specific suite
  call test_suite_init('Achilles Fortran Test Suite', test_suite)
  ! Test spectral function interface
  call test_spectral(test_suite)

  ! report the complete suite
  call test_suite_report(test_suite)
  ! finalize
  call test_suite_final(test_suite)

end program fortran_test
