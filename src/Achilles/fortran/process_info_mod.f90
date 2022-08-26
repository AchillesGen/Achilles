module libprocess_info
    use iso_c_binding

    private
    public :: process_info

    include "process_info_cdef.f90"

    type process_info 
        private
        type(c_ptr), pointer :: ptr
    contains
        procedure :: self => get_ptr
        procedure :: ids => process_ids 
        procedure :: multiplicity => process_multiplicity 
        procedure :: masses => process_masses 
        procedure :: add_state => process_add_state
    end type process_info

    interface process_info
        module procedure copy_process_info
    end interface

contains

    function copy_process_info(info)
        implicit none
        type(c_ptr), intent(in), target :: info
        type(process_info) :: copy_process_info
        copy_process_info%ptr => info
    end function

    function get_ptr(this)
        implicit none
        class(process_info), intent(in) :: this
        type(c_ptr) :: get_ptr
        get_ptr = this%ptr
    end function

    subroutine process_ids(this, ids, len)
        implicit none
        class(process_info), intent(in) :: this
        integer(c_long), dimension(:) :: ids
        integer(c_size_t), value :: len

        call process_ids_c(this%ptr, ids, len)
    end subroutine

    function process_multiplicity(this)
        implicit none
        class(process_info), intent(in) :: this
        integer(c_size_t) :: process_multiplicity

        process_multiplicity = process_multiplicity_c(this%ptr)
    end function

    subroutine process_masses(this, masses, len)
        implicit none
        class(process_info), intent(in) :: this
        integer(c_size_t), value :: len
        real(c_double), dimension(:) :: masses 

        call process_masses_c(this%ptr, masses, len)
    end subroutine

    subroutine process_add_state(this, initial, final, in, out)
        implicit none
        class(process_info), intent(inout) :: this
        integer(c_size_t), intent(in) :: in, out
        integer(c_long), dimension(:), intent(in) :: initial(in), final(out)

        call process_add_state_c(this%ptr, initial, final, in, out)
    end subroutine

end module
