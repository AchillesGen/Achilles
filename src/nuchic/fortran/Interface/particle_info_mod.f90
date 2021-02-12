module libpartinfo
    use iso_c_binding

    private
    public :: pinfo

    include "particle_info_cdef.f90"

    type pinfo
        private
        type(c_ptr) :: ptr ! Pointer to particle info obj
    contains
        ! Bind some functions to the type for cleaner syntax
        final :: delete_pinfo

        ! Member functions
        procedure :: self => get_ptr
        procedure :: name => get_name
        procedure :: pid => get_pid
        procedure :: charge => get_charge
        procedure :: spin => get_spin
        procedure :: mass => get_mass
        procedure :: width => get_width
    end type pinfo

    ! This will act as the constructor
    interface pinfo
        module procedure create_pinfo
        module procedure copy_constructor
    end interface

contains
    function create_pinfo(id)
        implicit none
        type(pinfo) :: create_pinfo
        integer, intent(in) :: id
        create_pinfo%ptr = create_pinfo_c(id)
    end function

    function copy_constructor(other)
        implicit none
        type(c_ptr), intent(in) :: other
        type(pinfo) :: copy_constructor
        copy_constructor%ptr = other
    end function

    subroutine delete_pinfo(this)
        implicit none
        type(pinfo) :: this
        call delete_pinfo_c(this%ptr)
    end subroutine

    function get_ptr(this)
        implicit none
        class(pinfo), intent(in) :: this
        type(c_ptr) :: get_ptr
        get_ptr = this%ptr
    end function

    function get_name(this)
        use libutilities
        implicit none
        character(len=:), allocatable :: get_name
        class(pinfo) :: this
        get_name = c2fstring(name_c(this%ptr))
    end function

    function get_pid(this)
        implicit none
        integer :: get_pid
        class(pinfo) :: this
        get_pid = pid_c(this%ptr)
    end function

    function get_charge(this)
        implicit none
        double precision :: get_charge
        class(pinfo) :: this
        get_charge = charge_c(this%ptr)
    end function

    function get_spin(this)
        implicit none
        double precision :: get_spin
        class(pinfo) :: this
        get_spin = spin_c(this%ptr)
    end function

    function get_mass(this)
        implicit none
        double precision :: get_mass
        class(pinfo) :: this
        get_mass = mass_c(this%ptr)
    end function

    function get_width(this)
        implicit none
        double precision :: get_width
        class(pinfo) :: this
        get_width = width_c(this%ptr)
    end function
end module
