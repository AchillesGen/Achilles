module libspectral_function
    use iso_c_binding

    private 
    public :: spectral_function

    include "spectral_function_cdef.f90"

    type spectral_function
        private
        type(c_ptr) :: ptr
    contains
        final :: delete_spectral

        procedure :: normalization => spectral_normalization
        procedure :: call => spectral_call
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
    end function

    subroutine delete_spectral(self)
        implicit none
        type(spectral_function) :: self
        call delete_spectral_function_c(self%ptr)
    end subroutine

    function spectral_normalization(self)
        implicit none
        class(spectral_function), intent(inout) :: self
        double precision :: spectral_normalization
        spectral_normalization = spectral_normalization_c(self%ptr)
    end function

    function spectral_call(self, p, e)
        implicit none
        class(spectral_function), intent(inout) :: self
        double precision, intent(in), value :: p, e
        double precision :: spectral_call
        spectral_call = spectral_call_c(self%ptr, p, e)
    end function
end module
