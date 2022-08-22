module libnucleus
    use iso_c_binding

    private 
    public :: nucleus

    include "nucleus_cdef.f90"

    type nucleus 
        private
        type(c_ptr) :: ptr
    contains
        procedure :: self => get_ptr
        procedure :: nucleon => get_nucleon 
        procedure :: add => add_particle 
    end type nucleus

contains

    function get_ptr(this)
        implicit none
        class(nucleus), intent(in) :: this
        type(c_ptr) :: get_ptr
        get_ptr = this%ptr
    end function

    function get_nucleon(this, i)
        use libparticle
        implicit none
        class(nucleus), intent(in) :: this
        integer(c_size_t), intent(in) :: i
        type(particle) :: get_nucleon

        call get_nucleon_c(this%ptr, i, get_nucleon%self())
    end function

    subroutine add_particle(this, part)
        use libparticle 
        implicit none
        class(nucleus), intent(in) :: this
        class(particle), intent(in) :: part 

        call add_particle_c(this%ptr, part%self())
    end subroutine
end module
