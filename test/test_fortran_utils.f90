module test_fortran_utils
    use libutilities
    use iso_c_binding

    private
    public :: test_string_conversion

contains
    
    subroutine test_string_conversion(test_suite)
        use unit_test
        type(test_suite_type) :: test_suite
        type(c_ptr) :: cstr
        character(len=:), allocatable :: fstr
        character(len=100) :: test_input

        call test_case_create('String conversion', test_suite)
        test_input = "test_filename.txt"
        cstr = f2cstring(trim(test_input))
        fstr = c2fstring(cstr)

        call assert_true(trim(test_input) == trim(fstr), suite=test_suite)
    end subroutine
end module
