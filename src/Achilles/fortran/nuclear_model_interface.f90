module nuclear_model_interface
    use nuclear_model_factory
    use nuclear_model
    use qe_spectral_model_test
    use res_spectral_model_test

    implicit none

    public
    private :: model_ptr

    class(model), pointer :: model_ptr
    type(ModelFactory) :: factory

contains

    subroutine register_all() bind(C, name="RegisterAll")
        ! Add all nuclear models to be registered here
        ! along with the function needed to build the model
        call factory%register_model("QE_Spectral_Func", build_qe_spec)
        call factory%register_model("RES_Spectral_Func", build_res_spec)
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

    function get_ps_name() bind(C, name="GetName_")
        use iso_c_binding
        use libutilities
        implicit none
        type(c_ptr) :: get_ps_name
        type(c_ptr) :: tmp
        character(len=:), allocatable :: fname

        fname = model_ptr%ps_name()
        get_ps_name = f2cstring(fname)
    end function

    subroutine currents(pids_in, pids_out, moms, nin, nout, qvec, ff, cur, nspin, nlorentz) bind(C, name="GetCurrents")
        use iso_c_binding
        use libvectors
        use libmap
        implicit none

        real(c_double), intent(in), dimension(4), target :: qvec
        integer(c_size_t), intent(in), value :: nin, nout, nspin, nlorentz
        complex(c_double_complex), dimension(nlorentz, nspin), intent(out) :: cur
        type(c_ptr), intent(in), value, target :: ff
        ! C++ is row major, while Fortran is column major
        ! This means the moms passed in are transposed
        integer(c_long), dimension(nin), intent(in) :: pids_in
        integer(c_long), dimension(nout), intent(in) :: pids_out
        real(c_double), intent(in), dimension(4, nin+nout) :: moms
        type(fourvector), dimension(nin) :: mom_in
        type(fourvector), dimension(nout) :: mom_out
        type(fourvector) :: qvector
        type(complex_map) :: ff_dict
        integer(c_size_t) :: i

        do i=1,nin
            mom_in(i) = fourvector(moms(1, i), moms(2, i), moms(3, i), moms(4, i))
        enddo

        do i=1,nout
            mom_out(i) = fourvector(moms(1, nin+i), moms(2, nin+i), moms(3, nin+i), moms(4, nin+i))
        enddo

        qvector = fourvector(qvec(1), qvec(2), qvec(3), qvec(4))
        ff_dict = complex_map(ff)
        call model_ptr%currents(pids_in, mom_in, nin, pids_out, mom_out, nout, qvector, ff_dict, cur, nspin, nlorentz)
    end subroutine

    function get_init_wgt(pids, pin, nin, nproton, nneutron) result(wgt) bind(c, name="GetInitialStateWeight")
        use iso_c_binding
        use libvectors
        implicit none

        integer(c_long), dimension(nin), intent(in) :: pids
        real(c_double), dimension(4, nin), intent(in) :: pin
        type(fourvector), dimension(nin) :: mom
        integer(c_size_t), intent(in), value :: nin, nproton, nneutron
        real(c_double) :: wgt
        integer(c_size_t) :: i

        do i=1,nin
            mom(i) = fourvector(pin(1, i), pin(2, i), pin(3, i), pin(4, i))
        enddo

        wgt = model_ptr%init_wgt(pids, mom, nin, nproton, nneutron)
    end function get_init_wgt

end module
