interface
    
    function create_interp1d_c(x, y, n) bind(C, name="CreateInterp1D")
        use iso_c_binding
        implicit none
        type(c_ptr) :: create_interp1d_c
        integer(c_int), intent(in), value :: n 
        real(c_double), dimension(n), intent(in) :: x
        real(c_double), dimension(n), intent(in) :: y
    end function

    subroutine delete_interp1d_c(self) bind(C, name="DeleteInterp1d")
        use iso_c_binding
        implicit none
        type(c_ptr) :: self
    end subroutine

    function interp1d_min_c(self) bind(C, name="Interp1DMin")
        use iso_c_binding
        implicit none
        type(c_ptr), intent(in), value :: self
        real(c_double) :: interp1d_min_c
    end function

    function interp1d_max_c(self) bind(C, name="Interp1DMax")
        use iso_c_binding
        implicit none
        type(c_ptr), intent(in), value :: self
        real(c_double) :: interp1d_max_c
    end function

    function interpolate1d_c(self, x) bind(C, name="Interpolate1D")
        use iso_c_binding
        implicit none
        type(c_ptr), intent(in), value :: self
        real(c_double), intent(in), value :: x
        real(c_double) :: interpolate1d_c
    end function

    function create_interp2d_c(x, y, z, n1, n2) bind(C, name="CreateInterp2D")
        use iso_c_binding
        implicit none
        type(c_ptr) :: create_interp2d_c
        integer(c_int), intent(in), value :: n1, n2
        real(c_double), dimension(n1), intent(in) :: x
        real(c_double), dimension(n2), intent(in) :: y
        real(c_double), dimension(n1*n2), intent(in) :: z
    end function

    subroutine delete_interp2d_c(self) bind(C, name="DeleteInterp2D")
        use iso_c_binding
        implicit none
        type(c_ptr), intent(in), value :: self
    end subroutine

    function interp2d_xmin_c(self) bind(C, name="Interp2DXMin")
        use iso_c_binding
        implicit none
        type(c_ptr), intent(in), value :: self
        real(c_double) :: interp2d_xmin_c
    end function

    function interp2d_xmax_c(self) bind(C, name="Interp2DXMax")
        use iso_c_binding
        implicit none
        type(c_ptr), intent(in), value :: self
        real(c_double) :: interp2d_xmax_c
    end function

    function interp2d_ymin_c(self) bind(C, name="Interp2DYMin")
        use iso_c_binding
        implicit none
        type(c_ptr), intent(in), value :: self
        real(c_double) :: interp2d_ymin_c
    end function

    function interp2d_ymax_c(self) bind(C, name="Interp2DYMax")
        use iso_c_binding
        implicit none
        type(c_ptr), intent(in), value :: self
        real(c_double) :: interp2d_ymax_c
    end function

    function interpolate2d_c(self, x, y) bind(C, name="Interpolate2D")
        use iso_c_binding
        implicit none
        type(c_ptr), intent(in), value :: self
        real(c_double), intent(in), value :: x, y
        real(c_double) :: interpolate2d_c
    end function

end interface
