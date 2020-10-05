module libinteraction
    use iso_c_binding

    private
    public :: interaction

    include "interaction_cdef.f90"

    type IFactory
        private
            character(len=20) :: factory_type
            class(interaction), pointer :: interaction_type
        contains
            procedure :: init
            procedure :: create_interaction
            procedure :: final
    end type IFactory

    type, abstract :: interaction
        contains
            procedure(gen_cross_section), deferred :: cross_section
            procedure(gen_make_momentum), deferred :: make_momentum
    end type interaction

    abstract interface
        function gen_cross_section(self, part1, part2)
            use iso_c_binding
            use libparticle
            import :: interaction
            implicit none
            class(interaction), intent(in) :: self
            type(particle), intent(in), value :: part1, part2
            real(c_double) :: gen_cross_section
        end function

        function gen_make_momentum(self, sameid, pcm, rans)
            use iso_c_binding
            use libvectors
            import :: interaction
            implicit none
            class(interaction), intent(in) :: self
            logical, intent(in) :: sameid
            real(c_double), intent(in) :: pcm
            real(c_double), dimension(2), intent(in) :: rans
            type(threevector) :: gen_make_momentum
        end function
    end interface

    type, extends(interaction) :: ConstantInteraction
        contains
            procedure :: cross_section => constant_xsec
            procedure :: make_momentum => constant_momentum
    end type ConstantInteraction

    class(interaction), pointer :: model

contains

    subroutine init(self, string)
        class(IFactory), intent(inout) :: self
        character(len=*), intent(in) :: string
        self%factory_type = trim(string)
        self%interaction_type => null()
    end subroutine init

    subroutine final(self)
        class(IFactory), intent(inout) :: self
        deallocate(self%interaction_type)
        nullify(self%interaction_type)
    end subroutine final

    function create_interaction(self) result(ptr)
        class(IFactory) :: self
        class(interaction), pointer :: ptr

        if(self%factory_type == "Constant") then
            if(associated(self%interaction_type)) deallocate(self%interaction_type)
            allocate(ConstantInteraction :: self%interaction_type)
            ptr => self%interaction_type
        end if
    end function create_interaction

    function constant_xsec(self, part1, part2)
        use libparticle
        class(ConstantInteraction), intent(in) :: self
        type(particle), intent(in), value :: part1, part2
        real(c_double) :: constant_xsec

        constant_xsec = 10
    end function

    function constant_momentum(self, sameid, pcm, rans)
        use iso_c_binding
        use libvectors
        class(ConstantInteraction), intent(in) :: self
        logical, intent(in) :: sameid
        real(c_double), intent(in) :: pcm
        real(c_double), dimension(2), intent(in) :: rans
        real(c_double) :: ctheta, stheta, phi
        type(threevector) :: constant_momentum

        ctheta = 2*rans(1)-1
        stheta = dsqrt(1-ctheta*ctheta)
        phi = 2*pi*rans(2)
        constant_momentum = threevector(pcm*stheta*dcos(phi), pcm*stheta*dsin(phi), pcm*ctheta)
    end function

    subroutine init_interaction(name) bind(C, name="InitializeInteraction")
        use iso_c_binding
        use libutilities
        implicit none
        type(c_ptr), intent(in), value :: name
        character(len=:), allocatable :: fname
        type(IFactory) :: factory

        fname = c2fstring(name)

        call factory%init(fname)
        model => factory%create_interaction()
    end subroutine

    function cross_section(part1, part2) bind(C, name="CrossSectionFortran")
        use iso_c_binding
        use libparticle
        use liblogging
        implicit none
        type(c_ptr), intent(in), value :: part1, part2
        type(particle) :: fpart1, fpart2
        real(c_double) :: cross_section

        fpart1 = particle(part1)
        fpart2 = particle(part2)

        cross_section = model%cross_section(fpart1, fpart2)
    end function

    function make_momentum(sameid, pcm, rans) bind(C, name="MakeMomentumFortran")
        use iso_c_binding
        use liblogging
        use libvectors
        implicit none
        logical, intent(in), value :: sameid
        real(c_double), intent(in), value :: pcm
        real(c_double), intent(in), dimension(2) :: rans
        type(threevector) :: mom
        type(c_ptr) :: make_momentum

        mom = model%make_momentum(sameid, pcm, rans)
        make_momentum = mom%copy()
    end function

end module
