module test_spectral_interface
    use libspectral_function

    private
    public :: test_spectral
    type(spectral_function) :: spectral

contains

    subroutine test_spectral(test_suite)
        use unit_test
        type(test_suite_type) :: test_suite
        double precision :: norm, val

        call test_case_create('Fortran Spectral function interface', test_suite)
        call assert_true(test_spectral_init("data/pke12_tot.data"), __FILE__, __LINE__, test_suite)
        norm = test_spectral_normalization()
        call assert_approximate(norm, 5.999988260179247d0, suite=test_suite)
        val = test_spectral_call(10d0, 22.5d0)
        call assert_approximate(val, 0.254d-8, suite=test_suite)
    end subroutine

    function test_spectral_init(filename) 
        use iso_c_binding
        implicit none
        character(len=*), intent(in) :: filename
        logical :: test_spectral_init

        spectral = spectral_function(filename)
        test_spectral_init = c_associated(spectral%self())
    end function

    function test_spectral_normalization()
        use iso_c_binding
        implicit none
        double precision :: test_spectral_normalization

        test_spectral_normalization = spectral%normalization()
    end function

    function test_spectral_call(p, E)
        use iso_c_binding
        implicit none
        double precision :: test_spectral_call
        double precision, intent(in) :: p, E

        test_spectral_call = spectral%call(p, E)
    end function

end module
