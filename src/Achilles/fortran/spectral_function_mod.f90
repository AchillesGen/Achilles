module libspectral_function
    use iso_c_binding

    private 
    public :: spectral_function

    include "spectral_function_cdef.f90"

    type spectral_function
        private
        type(c_ptr) :: ptr
    contains
        procedure :: normalization => spectral_normalization
        procedure :: callspectral => spectral_call
        procedure :: callmom => mom_density_call
        generic :: call => callspectral, callmom
        procedure :: self => spectral_self
    end type

    interface spectral_function
        module procedure create_spectral
    end interface

contains

    function create_spectral(filename)
        use libutilities
        implicit none
        character(len=*), intent(in) :: filename
        type(c_ptr) :: cfilename
        type(spectral_function) :: create_spectral

        cfilename = f2cstring(filename)
        create_spectral%ptr = create_spectral_function_c(cfilename)
        call free(cfilename)
    end function

    subroutine delete_spectral(self)
        implicit none
        type(spectral_function) :: self
        call delete_spectral_function_c(self%ptr)
    end subroutine

    function spectral_self(self)
        implicit none
        class(spectral_function), intent(inout) :: self
        type(c_ptr) :: spectral_self

        spectral_self = self%ptr
    end function

    function spectral_normalization(self)
        implicit none
        class(spectral_function), intent(in) :: self
        double precision :: spectral_normalization
        spectral_normalization = spectral_normalization_c(self%ptr)
    end function

    function spectral_call(self, p, e)
        implicit none
        class(spectral_function), intent(in) :: self
        double precision, intent(in), value :: p, e
        double precision :: spectral_call
        spectral_call = spectral_call_c(self%ptr, p, e)
    end function


    function mom_density_call(self, p)
        implicit none
        class(spectral_function), intent(in) :: self
        double precision, intent(in), value :: p
        double precision :: mom_density_call
        mom_density_call = mom_density_call_c(self%ptr, p)
    end function
end module
