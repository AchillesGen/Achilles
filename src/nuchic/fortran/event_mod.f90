module libevent
    use iso_c_binding

    private 
    public :: event

    include "event_cdef.f90"

    type event 
        private 
        type(c_ptr) :: ptr ! Pointer to the event object
    contains
        ! Bind some functions to the type for cleaner syntax
        procedure :: self => get_ptr
        procedure :: momentum => event_momentum
        procedure :: matrix_wgt => event_matrix_wgt 
        procedure :: total_xsec => event_total_cross_section 
        procedure :: select_nucleon => event_select_nucleon 
        procedure :: current_nucleus => event_nucleus 
        procedure :: set_mewgt => set_event_mewgt 
    end type event

    interface event
        module procedure copy_event
    end interface

contains

    function copy_event(evt)
        implicit none
        type(c_ptr), intent(in) :: evt
        type(event) :: copy_event 
        copy_event%ptr = evt
    end function

    function get_ptr(this)
        implicit none
        class(event), intent(in) :: this
        type(c_ptr) :: get_ptr
        get_ptr = this%ptr 
    end function

    function event_momentum(this, i)
        use libvectors 
        implicit none

        class(event), intent(in) :: this
        type(fourvector) :: event_momentum 
        integer(c_size_t), intent(in) :: i

        call event_momentum_c(this%ptr, i, event_momentum%self())
    end function
    
    subroutine event_matrix_wgt(this, i, wgt)
        implicit none

        class(event), intent(in) :: this
        integer(c_size_t), intent(in) :: i
        real(c_double), intent(in) :: wgt

        call event_matrix_wgt_c(this%ptr, i, wgt)
    end subroutine

    function event_total_cross_section(this)
        implicit none

        class(event), intent(in) :: this
        double precision :: event_total_cross_section 

        event_total_cross_section = event_total_cross_section_c(this%ptr)
    end function

    function event_select_nucleon(this)
        implicit none

        class(event), intent(in) :: this
        integer :: event_select_nucleon

        event_select_nucleon = event_select_nucleon_c(this%ptr)
    end function

    function event_nucleus(this)
        use libnucleus
        implicit none

        class(event), intent(in) :: this
        type(nucleus) :: event_nucleus

        call event_nucleus_c(this%ptr, event_nucleus%self())
    end function

    subroutine set_event_mewgt(this, wgt)
        implicit none

        class(event), intent(in) :: this
        double precision, intent(in) :: wgt

        call set_event_mewgt_c(this%ptr, wgt)
    end subroutine
end module
