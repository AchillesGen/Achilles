module factory
    public :: interaction, Create

    private

    type, abstract :: interaction
        character(100) :: value
    contains
        procedure(gen_xsec), deferred :: CrossSection
    end type

    abstract interface
        subroutine gen_xsec(this)
            import interaction
            class(interaction), intent(in) :: this
        end subroutine
    end interface

    type, extends(interaction) :: ConstantInteractions
    contains
        procedure :: CrossSection => ConstantCrossSection
    end type

contains
    subroutine ConstantCrossSection(this)
        class(ConstantInteractions), intent(in) :: this
        print *, "This is a " // this%value
    end subroutine

    subroutine Create(v, e)
        character(*), intent(in) :: v
        class(interaction), allocatable, intent(out) :: e

        if(v .EQ. 'ConstantInteractions') then
            allocate(e, source=ConstantInteractions('ConstantInteractions'))
        endif
    end subroutine Create
end module factory
