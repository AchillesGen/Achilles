program fortran_test
  use unit_test
  use test_spectral_interface

  implicit none

  type(test_suite_type) :: test_suite
  double precision :: norm, val

  ! example with specific suite
  call test_suite_init('Achilles Fortran Test Suite', test_suite)
  call test_case_create('Fortran Spectral function interface', test_suite)

  call assert_true(test_spectral_init("data/pke12_tot.data"), __FILE__, __LINE__, test_suite)
  norm = test_spectral_normalization()
  call assert_approximate(norm, 5.999988260179247d0, suite=test_suite)
  val = test_spectral_call(10d0, 22.5d0)
  call assert_approximate(val, 0.254d-8, suite=test_suite)

  ! report the complete suite
  call test_suite_report(test_suite)
  ! finalize
  call test_suite_final(test_suite)

end program fortran_test
