module libparticle
    use iso_c_binding

    private
    public :: particle

    include "particle_cdef.f90"

    type particle
        private
        type(c_ptr) :: ptr
    contains
        final :: delete_particle

        ! Member functions
        procedure :: status => get_status
        procedure :: info => get_info
        procedure :: momentum => get_momentum
        procedure :: position => get_position
        procedure :: set_momentum => set_particle_momentum
        procedure :: set_position => set_particle_position
    end type particle

    interface particle
        module procedure create_particle
    end interface

contains
    function create_particle(id, momentum, position, status)
        use libvectors
        implicit none
        type(particle) :: create_particle
        integer, intent(in) :: id
        type(fourvector), intent(in) :: momentum
        type(threevector), intent(in) :: position
        integer, intent(in) :: status
        create_particle%ptr = create_particle_c(id, momentum%self(), position%self(), status)
    end function

    subroutine delete_particle(this)
        implicit none
        type(particle) :: this
        call delete_particle_c(this%ptr)
    end subroutine

    function get_status(this)
        implicit none
        integer :: get_status
        class(particle) :: this
        get_status = particle_status_c(this%ptr)
    end function

    function get_info(this)
        use libpartinfo
        implicit none
        type(pinfo) :: get_info
        class(particle) :: this
        get_info = pinfo(particle_info_c(this%ptr))
    end function

    function get_momentum(this)
        use libvectors
        implicit none
        type(fourvector) :: get_momentum
        class(particle) :: this
        get_momentum = fourvector(particle_momentum_c(this%ptr))
    end function

    function get_position(this)
        use libvectors
        implicit none
        type(threevector) :: get_position
        class(particle) :: this
        get_position = threevector(particle_position_c(this%ptr))
    end function

    subroutine set_particle_momentum(this, momentum)
        use libvectors
        implicit none
        type(fourvector), intent(in) :: momentum
        class(particle) :: this
        this%ptr = set_particle_momentum_c(this%ptr, momentum%self())
    end subroutine

    subroutine set_particle_position(this, position)
        use libvectors
        implicit none
        type(threevector), intent(in) :: position
        class(particle), intent(inout) :: this
        this%ptr = set_particle_position_c(this%ptr, position%self())
    end subroutine

end module
