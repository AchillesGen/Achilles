program fortran_test
  use unit_test

  implicit none

  type(test_suite_type) :: test_suite

  ! example with specific suite
  call test_suite_init('Nuchic Fortran Test Suite', test_suite)
  call test_case_create('Dummy Test', test_suite)
  ! suite = SUITE need in this case (cause optional argument eps, file_name, line_number is missing)
  call assert_approximate(1.0, 1.0, __FILE__, __LINE__, suite=test_suite)

  call test_case_create('Fortran YAML', test_suite)

  ! report the complete suite
  call test_suite_report(test_suite)
  ! finalize
  call test_suite_final(test_suite)

end program fortran_test
