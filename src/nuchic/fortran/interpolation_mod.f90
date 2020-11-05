module libinterpolate
    use iso_c_binding

    private
    public :: interp1d, interp2d

    include "interpolation_cdef.f90"

    type interp1d
        private
        type(c_ptr) :: ptr
    contains
        final :: delete_interp1d

        procedure :: min => interp1d_min
        procedure :: max => interp1d_max
        procedure :: call => interp1d_call
    end type

    interface interp1d
        module procedure create_interp1d
    end interface

    type interp2d
        private
        type(c_ptr) :: ptr
    contains
        final :: delete_interp2d

        procedure :: xmin => interp2d_xmin
        procedure :: xmax => interp2d_xmax
        procedure :: ymin => interp2d_ymin
        procedure :: ymax => interp2d_ymax
        procedure :: call => interp2d_call
    end type

    interface interp2d
        module procedure create_interp2d
        module procedure create_interp2d2
    end interface

contains

    function create_interp1d(x, y, n, mode)
        implicit none
        integer, intent(in) :: n, mode
        double precision, dimension(n), intent(in) :: x, y
        type(interp1d) :: create_interp1d
        create_interp1d%ptr = create_interp1d_c(x, y, n, mode)
    end function

    subroutine delete_interp1d(this)
        implicit none
        type(interp1d) :: this
        call delete_interp1d_c(this%ptr)
    end subroutine

    function interp1d_min(this)
        implicit none
        class(interp1d), intent(in) :: this
        double precision :: interp1d_min
        interp1d_min = interp1d_min_c(this%ptr)
    end function

    function interp1d_max(this)
        implicit none
        class(interp1d), intent(in) :: this
        double precision :: interp1d_max
        interp1d_max = interp1d_max_c(this%ptr)
    end function

    function interp1d_call(this, x)
        implicit none
        class(interp1d), intent(in) :: this
        double precision, intent(in), value :: x
        double precision :: interp1d_call
        interp1d_call = interpolate1d_c(this%ptr, x)
    end function

    function create_interp2d(x, y, z, n1, n2, mode)
        implicit none
        integer, intent(in) :: n1, n2, mode
        double precision, dimension(n1), intent(in) :: x
        double precision, dimension(n2), intent(in) :: y
        double precision, dimension(n1*n2), intent(in) :: z
        type(interp2d) :: create_interp2d
        create_interp2d%ptr = create_interp2d_c(x, y, z, n1, n2, mode)
    end function

    function create_interp2d2(x, y, z, n1, n2, mode)
        implicit none
        integer, intent(in) :: n1, n2, mode
        double precision, dimension(n1), intent(in) :: x
        double precision, dimension(n2), intent(in) :: y
        double precision, dimension(n1,n2), intent(in) :: z
        type(interp2d) :: create_interp2d2
        create_interp2d2%ptr = create_interp2d_c(x, y, transpose(z), n1, n2, mode)
    end function

    subroutine delete_interp2d(this)
        implicit none
        type(interp2d) :: this
        call delete_interp1d_c(this%ptr)
    end subroutine

    function interp2d_xmin(this)
        implicit none
        class(interp2d), intent(in) :: this
        double precision :: interp2d_xmin
        interp2d_xmin = interp2d_xmin_c(this%ptr)
    end function

    function interp2d_xmax(this)
        implicit none
        class(interp2d), intent(in) :: this
        double precision :: interp2d_xmax
        interp2d_xmax = interp2d_xmax_c(this%ptr)
    end function

    function interp2d_ymin(this)
        implicit none
        class(interp2d), intent(in) :: this
        double precision :: interp2d_ymin
        interp2d_ymin = interp2d_ymin_c(this%ptr)
    end function

    function interp2d_ymax(this)
        implicit none
        class(interp2d), intent(in) :: this
        double precision :: interp2d_ymax
        interp2d_ymax = interp2d_ymax_c(this%ptr)
    end function

    function interp2d_call(this, x, y)
        implicit none
        class(interp2d), intent(in) :: this
        double precision, intent(in), value :: x, y
        double precision :: interp2d_call
        interp2d_call = interpolate2d_c(this%ptr, x, y)
    end function

end module
