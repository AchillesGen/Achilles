module nuclear_model_factory
    use nuclear_model
    implicit none

    private 
    public :: ModelFactory

    abstract interface
        function constructor()
            import model 
            class(model), pointer :: constructor
        end function constructor
    end interface

    type :: ModelConstructor
        character(len=20) :: name
        procedure(constructor), pointer, nopass :: construct 
    end type ModelConstructor

    type :: ModelFactory
            class(model), pointer :: model_ptr
            type(ModelConstructor), allocatable :: model_list(:)
        contains
            procedure :: init
            procedure :: create_model
            procedure :: print_models
            procedure :: register_model
            procedure :: final
    end type ModelFactory

contains

    subroutine init(self)
        class(ModelFactory), intent(inout) :: self
        self % model_ptr => null()
    end subroutine init

    subroutine final(self)
        class(ModelFactory), intent(inout) :: self
        deallocate(self % model_ptr)
        nullify(self % model_ptr)
    end subroutine final

    function create_model(self, name) result(ptr)
        class(ModelFactory) :: self
        character(len=*), intent(in) :: name
        class(model), pointer :: ptr
        integer :: i

        if(allocated(self % model_list)) then
            do i = 1, size(self % model_list)
                if(trim(name) == trim(self % model_list(i) % name)) then
                    self % model_ptr => self % model_list(i) % construct()
                    ptr => self % model_ptr
                    return 
                endif
            enddo
        endif

        ! Pass throw to Achilles C++
        ! write(*,*) "FortranModel: Model ", trim(name), " is undefined"
        ! stop 
    end function create_model

    subroutine print_models(self)
        class(ModelFactory), intent(inout) :: self
        integer :: i

        if(allocated(self % model_list)) then
            do i = 1, size(self % model_list)
                write(*,*) " - ", trim(self % model_list(i) % name)
            enddo
        endif
    end subroutine print_models

    subroutine register_model(self, name, construct)
        class(ModelFactory), intent(inout) :: self
        character(len=*), intent(in) :: name
        procedure(constructor) :: construct
        type(ModelConstructor) :: build_model
        integer :: i

        ! Ensure name is not a duplicate
        if(allocated(self % model_list)) then
            do i = 1, size(self % model_list)
                if(name == self % model_list(i) % name) then
                    ! Pass throw to Achilles C++
                    write(*, *) "FortranModel: Model named ", trim(name), " already exists!"
                    stop 
                endif
            enddo
        endif

        build_model % construct => construct 
        build_model % name = trim(name)

        if(.not. allocated(self % model_list)) allocate(self % model_list(0))
        self % model_list = [self % model_list, build_model]
    end subroutine
end module nuclear_model_factory
