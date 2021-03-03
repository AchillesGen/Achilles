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
            procedure(gen_cross_sections), deferred :: cross_sections
            procedure(gen_final_state), deferred :: final_state
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

        subroutine gen_cross_sections(self, part1, part2, results, length)
            use iso_c_binding
            use libparticle
            import :: interaction
            implicit none
            class(interaction), intent(in) :: self
            type(particle), intent(in), value :: part1, part2
            integer(c_int), intent(out) :: length
            real(c_double), intent(out), dimension(:), pointer :: results
        end subroutine

        subroutine gen_final_state(self, part1, part2, results, length)
            use iso_c_binding
            use libparticle
            import :: interaction
            implicit none
            class(interaction), intent(in) :: self
            type(particle), intent(in), value :: part1, part2
            integer(c_int), intent(out) :: length
            type(particle), intent(out), dimension(:), pointer :: results
        end subroutine
    end interface

    type, extends(interaction) :: ConstantInteraction
        contains
            procedure :: cross_section => constant_xsec
            procedure :: cross_sections => constant_xsecs
            procedure :: final_state => constant_final_state
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

    subroutine constant_xsecs(self, part1, part2, results, length)
        use libparticle
        class(ConstantInteraction), intent(in) :: self
        type(particle), intent(in), value :: part1, part2
        integer(c_int), intent(out) :: length
        real(c_double), intent(out), dimension(:), pointer :: results

        length = 1
        allocate(results(length))

        results(1) = 10
    end subroutine

    subroutine constant_final_state(self, part1, part2, results, length)
        use libparticle
        class(ConstantInteraction), intent(in) :: self
        type(particle), intent(in), value :: part1, part2
        integer(c_int), intent(out) :: length
        type(particle), intent(out), dimension(:), pointer :: results

        length = 2
        allocate(results(length))
        results(1) = particle(part1%self())
        results(2) = particle(part2%self())
    end subroutine

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

    subroutine cross_sections(part1, part2, results, length) bind(C, name="CrossSectionsFortran")
        use iso_c_binding
        use libparticle
        implicit none
        type(c_ptr), intent(in), value :: part1, part2
        type(particle) :: fpart1, fpart2
        type(c_ptr), intent(out) :: results
        integer(c_int), intent(out) :: length
        double precision, dimension(:), pointer :: tmp_results

        fpart1 = particle(part1)
        fpart2 = particle(part2)

        call model%cross_sections(fpart1, fpart2, tmp_results, length)
        results = c_loc(tmp_results(1))
    end subroutine

    subroutine final_state(part1, part2, results, length) bind(C, name="GenerateFinalStateFortran")
        use iso_c_binding
        use libparticle
        implicit none
        type(c_ptr), intent(in), value :: part1, part2
        type(particle) :: fpart1, fpart2
        type(c_ptr), intent(out) :: results
        integer(c_int), intent(out) :: length
        type(particle), dimension(:), pointer :: tmp_results

        fpart1 = particle(part1)
        fpart2 = particle(part2)

        call model%final_state(fpart1, fpart2, tmp_results, length)
        results = c_loc(tmp_results(1))
    end subroutine

end module
