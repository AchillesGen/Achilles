module libevent
    use iso_c_binding

    private
    public :: event

    type event
        private
        type(c_ptr) :: ptr ! Pointer to the event object
    contains
        ! Bind some functions to the type for cleaner syntax
        procedure :: self => get_ptr
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
end module
