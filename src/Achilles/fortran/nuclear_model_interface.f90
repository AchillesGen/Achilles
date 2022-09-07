module nuclear_model_interface
    use nuclear_model_factory
    use nuclear_model
    implicit none

    public
    private :: model_ptr

    class(model), pointer :: model_ptr
    type(ModelFactory) :: factory

contains

    subroutine register_all() bind(C, name="RegisterAll")
        ! Add all nuclear models to be registered here
        ! along with the function needed to build the model
        ! call factory%register_model("test", build_test)
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

    subroutine currents(pids_in, pids_out, moms, nin, nout, qvec, ff, len_ff, cur, nspin, nlorentz) bind(C, name="GetCurrents")
        use iso_c_binding
        use libvectors
        implicit none

        real(c_double), intent(in), dimension(4), target :: qvec
        integer(c_size_t), intent(in), value :: len_ff, nin, nout, nspin, nlorentz
        complex(c_double_complex), dimension(nlorentz, nspin), intent(out) :: cur
        complex(c_double_complex), dimension(len_ff), intent(in), target :: ff
        ! C++ is row major, while Fortran is column major
        ! This means the moms passed in are transposed
        integer(c_int), dimension(nin), intent(in) :: pids_in
        integer(c_int), dimension(nout), intent(in) :: pids_out
        real(c_double), intent(in), dimension(4, nin+nout) :: moms
        type(fourvector), dimension(nin) :: mom_in
        type(fourvector), dimension(nout) :: mom_out
        type(fourvector) :: qvector
        integer(c_size_t) :: i

        do i=1,nin
            mom_in(i) = fourvector(moms(1, i), moms(2, i), moms(3, i), moms(4, i))
        enddo

        do i=1,nout
            mom_out(i) = fourvector(moms(1, nin+i), moms(2, nin+i), moms(3, nin+i), moms(4, nin+i))
        enddo

        qvector = fourvector(qvec(1), qvec(2), qvec(3), qvec(4))

        call model_ptr%currents(pids_in, mom_in, nin, pids_out, mom_out, nout, qvector, ff, len_ff, cur, nspin, nlorentz)
    end subroutine

!    function allowed_states(cinfo) bind(C, name="GetAllowedStates")
!        use iso_c_binding
!        use libprocess_info
!        implicit none
!
!        type(c_ptr) :: cinfo
!        type(process_info) :: info
!        logical :: allowed_states
!
!        info = process_info(cinfo)
!        allowed_states = model_ptr%states(info) 
!    end function

!    function nspins() bind(C, name="GetNSpins")
!        use iso_c_binding
!        implicit none
!
!        integer(c_size_t) :: nspins
!
!        nspins = model_ptr%nspins()
!    end function

!    function fill_nucleus(evt, xsec, len) bind(C, name="FillNucleus")
!        use iso_c_binding
!        use libevent
!        implicit none
!
!        type(c_ptr) :: evt
!        type(event) :: fevt
!        integer(c_size_t), value :: len
!        real(c_double), dimension(len), intent(in) :: xsec
!        logical :: fill_nucleus
!
!        fevt = event(evt)
!        fill_nucleus = model_ptr%fill_nucleus(fevt, xsec, len)
!    end function
end module
