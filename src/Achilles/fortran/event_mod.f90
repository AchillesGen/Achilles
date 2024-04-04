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
        procedure :: current_nucleus => event_nucleus 
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

    function event_nucleus(this)
        use libnucleus
        implicit none

        class(event), intent(in) :: this
        type(nucleus) :: event_nucleus

        call event_nucleus_c(this%ptr, event_nucleus%self())
    end function
end module
