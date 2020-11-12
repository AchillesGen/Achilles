program good_test
  use unit_test

  implicit none

  type(test_suite_type) :: specific_suite

  ! example with specific suite
  call test_suite_init('my specific test suite', specific_suite)
  call test_case_create('Specific Test 1', specific_suite)
  ! suite = SUITE need in this case (cause optional argument eps, file_name, line_number is missing)
  call assert_approximate(1.0, 2.0, __FILE__, __LINE__, suite=specific_suite)

  ! report the complete suite
  call test_suite_report(specific_suite)
  ! finalize
  call test_suite_final(specific_suite)

end program good_test
