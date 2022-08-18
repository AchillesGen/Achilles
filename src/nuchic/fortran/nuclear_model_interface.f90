module nuclear_model_interface
    use nuclear_model_factory
    use nuclear_model
    implicit none

    public
    private :: model_ptr, factory

    class(model), pointer :: model_ptr
    type(ModelFactory) :: factory

    type, extends(model) :: test
        contains
            procedure :: init => test_init
            procedure :: fill_nucleus => test_fill
            procedure :: nspins => test_spins
            procedure :: states => test_states
            procedure :: currents => test_currents
            procedure :: ps_name => test_ps
            procedure :: mode => test_mode
            procedure :: cleanup => test_cleanup
    end type


contains

    subroutine register_all() bind(C, name="RegisterAll")
        ! Add all nuclear models to be registered here
        ! along with the function needed to build the model
        call factory%register_model("test", build_test)
        call factory%init()
    end subroutine

    subroutine print_all() bind(C, name="ListModels")
        call factory%print_models()
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
        integer(c_int), intent(in) :: size
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

    function get_ps_name() bind(C, name="GetName")
        use iso_c_binding
        use libutilities
        implicit none
        type(c_ptr) :: get_ps_name
        type(c_ptr) :: tmp
        character(len=:), allocatable :: fname

        fname = model_ptr%ps_name()
        get_ps_name = f2cstring(fname)
    end function

    subroutine currents(moms, nin, nout, qvec, ff, len_ff, cur, len_cur) bind(C, name="GetCurrents")
        use iso_c_binding
        use libvectors
        implicit none

        real(c_double), intent(in), dimension(4), target :: qvec
        integer(c_size_t), intent(in), value :: len_ff, nin, nout
        integer(c_int), intent(out) :: len_cur
        type(c_ptr), intent(out) :: cur
        complex(c_double_complex), dimension(len_ff), intent(in), target :: ff
        ! C++ is row major, while Fortran is column major
        ! This means the moms passed in are transposed
        real(c_double), intent(in), dimension(4, nin+nout) :: moms
        type(fourvector), dimension(nin) :: mom_in
        type(fourvector), dimension(nout) :: mom_out
        type(fourvector) :: qvector
        integer(c_size_t) :: i
        complex(c_double_complex), dimension(:), pointer :: tmp_results

        do i=1,nin
            mom_in(i) = fourvector(moms(1, i), moms(2, i), moms(3, i), moms(4, i))
        enddo

        do i=1,nout
            mom_out(i) = fourvector(moms(1, nin+i), moms(2, nin+i), moms(3, nin+i), moms(4, nin+i))
        enddo

        qvector = fourvector(qvec(1), qvec(2), qvec(3), qvec(4))

        call model_ptr%currents(mom_in, nin, mom_out, nout, qvector, ff, len_ff, tmp_results, len_cur)
        cur = c_loc(tmp_results(1))
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

    function test_init(self, filename)
        class(test), intent(inout) :: self
        character(len=*), intent(in) :: filename
        logical :: test_init
        write(*, *) filename
        test_init = .true.
    end function

    function build_test()
        class(model), pointer :: build_test
        allocate(test :: build_test)
    end function build_test

    subroutine test_cleanup(self)
        class(test), intent(inout) :: self
    end subroutine

    function test_mode(self)
        class(test), intent(inout) :: self
        integer :: test_mode
        test_mode = -1
    end function

    function test_ps(self)
        class(test), intent(inout) :: self
        character(len=:), allocatable :: test_ps
        test_ps = "QESpectral"
    end function

    subroutine test_currents(self, mom_in, nin, mom_out, nout, qvec, ff, len_ff, cur, len_cur)
        use iso_c_binding
        use libvectors
        class(test), intent(inout) :: self
        integer(c_size_t), intent(in), value :: nin, nout, len_ff
        integer(c_int), intent(out) :: len_cur
        complex(c_double_complex), dimension(len_ff), intent(in) :: ff
        type(fourvector) :: qvec
        type(fourvector), dimension(nin), intent(in) :: mom_in
        type(fourvector), dimension(nout), intent(in) :: mom_out
        complex(c_double_complex), dimension(:), intent(out), pointer :: cur

        len_cur = 8
        allocate(cur(len_cur))
        cur(1) = 1
        cur(2) = 2
        cur(3) = 3
        cur(4) = 4
        cur(5) = 5
        cur(6) = 6
        cur(7) = 7
        cur(8) = 8
    end subroutine

    function test_states(self, info)
        use libprocess_info
        use iso_c_binding
        class(test), intent(inout) :: self
        type(process_info), intent(inout) :: info
        logical :: test_states
        integer(c_size_t) :: in, out
        integer(c_long), dimension(:), allocatable :: initial, final
        in = 1
        out = 1
        allocate(initial(in))
        allocate(final(out))
        initial(1) = 2212
        final(1) = 2212
        test_states = .true.
        call info%add_state(initial, final, in, out)
        deallocate(initial)
        deallocate(final)
    end function

    function test_spins(self)
        use iso_c_binding
        class(test), intent(inout) :: self
        integer(c_size_t) test_spins
        test_spins = 2
    end function

    function test_fill(self, evt, xsec, len)
        use libevent
        use iso_c_binding
        class(test), intent(inout) :: self
        class(event), intent(inout) :: evt
        integer(c_size_t), value :: len
        real(c_double), dimension(len), intent(in) :: xsec
        logical :: test_fill
        test_fill = .true.
    end function
end module
