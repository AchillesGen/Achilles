interface

    function create_spectral_function_c(filename) bind(C, name="LoadSpectralFunction")
        use iso_c_binding
        implicit none
        type(c_ptr) :: create_spectral_function_c 
        type(c_ptr), intent(in) :: filename
    end function

    subroutine delete_spectral_function_c(self) bind(C, name="DeleteSpectralFunction")
        use iso_c_binding
        implicit none
        type(c_ptr), intent(inout) :: self
    end subroutine

    function spectral_normalization_c(self) bind(C, name="SpectralNormalization")
        use iso_c_binding
        implicit none
        type(c_ptr), intent(inout) :: self
        real(c_double) :: spectral_normalization_c
    end function

    function spectral_call_c(self, p, e) bind(C, name="SpectralFunction")
        use iso_c_binding
        implicit none
        type(c_ptr), intent(inout) :: self
        real(c_double), intent(in), value :: p, e
        real(c_double) :: spectral_call_c
    end function

end interface
