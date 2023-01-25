module test_spectral_interface
    use libspectral_function

    private
    public :: test_spectral_init, test_spectral_normalization, test_spectral_call
    type(spectral_function) :: spectral

contains

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
