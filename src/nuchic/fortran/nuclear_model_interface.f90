module nuclear_model_interface
    use nuclear_model_factory
    use nuclear_model
    implicit none

    public
    private :: model_ptr, factory

    class(model), pointer :: model_ptr
    type(ModelFactory) :: factory

contains

    subroutine register_all() bind(C, name="RegisterAll")
        ! Add all nuclear models to be registered here
        ! along with the function needed to build the model
        ! call factory%register_model("test", build_test)
        call factory%init()
    end subroutine

    function create_model(name) bind(C, name="CreateModel")
        use iso_c_binding
        use libutilities
        implicit none

        type(c_ptr), intent(in), value :: name
        logical :: create_model
        character(len=:), allocatable :: fname

        fname = c2fstring(name)

        model_ptr => factory%create_model(trim(fname))
      
        create_model = .false.
        if(associated(model_ptr)) create_model = .true.
    end function

    function init_model(name) bind(C, name="InitModel")
        use iso_c_binding
        use libutilities
        implicit none

        type(c_ptr), intent(in), value :: name
        logical :: init_model
        character(len=:), allocatable :: fname

        fname = c2fstring(name)

        init_model = model_ptr%init(fname)
    end function

    subroutine clean_event(cur, size) bind(C, name="CleanUpEvent")
        use iso_c_binding
        implicit none

        type(c_ptr), intent(in) :: cur 
        integer(c_size_t), intent(in) :: size
        complex(c_double_complex), dimension(:), pointer :: tmp_cur

        call c_f_pointer(cur, tmp_cur, [size])
        deallocate(tmp_cur)
        nullify(tmp_cur)
    end subroutine

    subroutine clean_model() bind(C, name="CleanUpModel")
        implicit none
        deallocate(model_ptr)
        nullify(model_ptr)
    end subroutine

    function get_mode() bind(C, name="GetMode")
        use iso_c_binding
        implicit none
        integer(c_int) :: get_mode

        get_mode = model_ptr%mode()
    end function

    subroutine get_ps_name(name) bind(C, name="GetName")
        use iso_c_binding
        use libutilities
        implicit none
        type(c_ptr), intent(out) :: name
        character(len=:), allocatable :: fname

        fname = model_ptr%ps_name()
        name = f2cstring(fname)
    end subroutine

    !TODO: Figure out the details of this function
    !TODO: Develop fortran form factor wrapper
    subroutine currents(evt, ff, len_ff, cur, len_cur) bind(C, name="GetCurrents")
        use iso_c_binding
        use libevent
        implicit none

        type(c_ptr) :: evt, ff
        integer(c_size_t) :: len_ff, len_cur
        complex(c_double_complex), dimension(:), allocatable :: cur
    end subroutine

    function allowed_states(cinfo) bind(C, name="GetAllowedStates")
        use iso_c_binding
        use libprocess_info
        implicit none

        type(c_ptr) :: cinfo
        type(process_info) :: info
        logical :: allowed_states

        info = process_info(cinfo)
        allowed_states = model_ptr%states(info) 
    end function

    function nspins() bind(C, name="GetNSpins")
        use iso_c_binding
        implicit none

        integer(c_size_t) :: nspins

        nspins = model_ptr%nspins()
    end function

    function fill_nucleus(evt, xsec, len) bind(C, name="FillNucleus")
        use iso_c_binding
        use libevent
        implicit none

        type(c_ptr) :: evt
        type(event) :: fevt
        integer(c_size_t), value :: len
        real(c_double), dimension(len), intent(in) :: xsec
        logical :: fill_nucleus

        fevt = event(evt)
        fill_nucleus = model_ptr%fill_nucleus(fevt, xsec, len)
    end function
end module
