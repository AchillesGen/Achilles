module nuclear_model_interface
    use iso_c_binding
    use nuclear_model_factory
    use nuclear_model
    use qe_spectral_model
    use res_spectral_model
    use intf_spectral_model

    implicit none

    public
    private :: model_holder, model_vect, models

    type(ModelFactory) :: factory

    type model_holder
        class(model), pointer :: model_ptr
    end type model_holder

    type model_vect
        private
        type(model_holder), allocatable  :: model_holder(:)
        integer(c_size_t) :: size, capacity
        logical :: initialized = .false.
        contains
            ! Member functions
            procedure :: init => model_vect_init
            procedure :: add_model => model_vect_add_model
            procedure :: get_model => model_vect_get_model
            procedure :: resize => model_vect_resize
    end type model_vect

    type(model_vect) :: models

contains

    subroutine model_vect_init(this)
        class(model_vect), intent(inout) :: this
        this%size = 0
        this%capacity = 1
        allocate(this%model_holder(this%capacity))
    end subroutine

    subroutine model_vect_add_model(this, model)
        class(model_vect), intent(inout) :: this
        type(model_holder), intent(in) :: model
        integer(c_size_t) :: i

        if(.not. this%initialized) then
            call this%init()
            this%initialized = .true.
        endif

        if(this%size == this%capacity) then
            call this%resize(2*this%capacity)
        endif

        i = this%size + 1
        this%model_holder(i) = model
        this%size = i
    end subroutine

    function model_vect_get_model(this, idx) result(model)
        class(model_vect), intent(in) :: this
        integer(c_size_t), intent(in) :: idx
        type(model_holder) :: model
        model = this%model_holder(idx)
    end function

    subroutine model_vect_resize(this, new_size)
        class(model_vect), intent(inout) :: this
        integer(c_size_t), intent(in), optional :: new_size
        integer :: i
        type(model_holder), allocatable :: tmp(:)

        if(allocated(this%model_holder)) then
            call move_alloc(this%model_holder, tmp)
        endif

        if(present(new_size)) then
            this%capacity = new_size
        else
            this%capacity = 2*this%capacity
        endif

        allocate(this%model_holder(this%capacity))

        if(allocated(tmp)) then
            do i=1,this%size
                this%model_holder(i) = tmp(i)
            enddo
        end if
    end subroutine

    subroutine register_all() bind(C, name="RegisterAll")
        type(qe_spec) :: qe
        type(res_spec) :: res
        type(intf_spec) :: intf
        ! Add all nuclear models to be registered here
        ! along with the function needed to build the model
        call factory%register_model(qe%model_name(), build_qe_spec)
        call factory%register_model(res%model_name(), build_res_spec)
        call factory%register_model(intf%model_name(), build_intf_spec)
        call factory%init()
    end subroutine

    subroutine print_all() bind(C, name="ListModels")
        call factory%print_models()
    end subroutine

    function create_model(name, idx) bind(C, name="CreateModel")
        use iso_c_binding
        use libutilities
        implicit none

        type(c_ptr), intent(in), value :: name
        integer(c_size_t), intent(inout) :: idx
        logical :: create_model
        character(len=:), allocatable :: fname
        class(model), pointer :: model_ptr
        type(model_holder) :: holder

        fname = c2fstring(name)

        model_ptr => factory%create_model(trim(fname))
      
        create_model = .false.
        if(associated(model_ptr)) then
            create_model = .true.
            holder%model_ptr => model_ptr
            call models%add_model(holder)
            idx = models%size
        endif
    end function

    function init_model(name, params, idx) bind(C, name="InitModel")
        use iso_c_binding
        use libutilities
        use libmap
        implicit none

        type(c_ptr), intent(in), value :: name
        integer(c_size_t), intent(in), value :: idx
        type(c_ptr), intent(in), value, target :: params
        type(c_ptr) :: model_cptr
        logical :: init_model
        character(len=:), allocatable :: fname
        type(map) :: params_dict
        type(model_holder) :: model

        model = models%get_model(idx)
        fname = c2fstring(name)
        params_dict = map(params)
        init_model = model%model_ptr%init(fname,params_dict)
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

    subroutine clean_model(idx) bind(C, name="CleanUpModel")
        implicit none
        integer(c_size_t), intent(in), value :: idx
        type(model_holder) :: model

        model = models%get_model(idx)
        deallocate(model%model_ptr)
        nullify(model%model_ptr)
    end subroutine

    function get_mode(idx) bind(C, name="GetMode")
        use iso_c_binding
        implicit none
        integer(c_size_t), intent(in), value :: idx
        integer(c_int) :: get_mode
        type(model_holder) :: model

        model = models%get_model(idx)
        get_mode = model%model_ptr%mode()
    end function

    function get_frame(idx) bind(C, name="GetFrame")
        use iso_c_binding
        implicit none
        integer(c_size_t), intent(in), value :: idx
        integer(c_int) :: get_frame
        type(model_holder) :: model

        model = models%get_model(idx)
        get_frame = model%model_ptr%frame()
    end function

    function get_model_name(idx) bind(C, name="ModelName")
        use iso_c_binding
        use libutilities
        implicit none
        integer(c_size_t), intent(in), value :: idx
        type(c_ptr) :: get_model_name
        type(c_ptr) :: tmp
        character(len=:), allocatable :: fname
        type(model_holder) :: model

        model = models%get_model(idx)
        fname = model%model_ptr%model_name()
        get_model_name = f2cstring(fname)
    end function

    function get_ps_name(idx) bind(C, name="GetName_")
        use iso_c_binding
        use libutilities
        implicit none
        integer(c_size_t), intent(in), value :: idx
        type(c_ptr) :: get_ps_name
        type(c_ptr) :: tmp
        character(len=:), allocatable :: fname
        type(model_holder) :: model

        model = models%get_model(idx)
        fname = model%model_ptr%ps_name()
        get_ps_name = f2cstring(fname)
    end function

    subroutine currents(idx, pids_in, pids_out, pids_spect, moms, nin, nout, nspect, qvec, ff, cur, nspin, nlorentz) bind(C, name="GetCurrents")
        use iso_c_binding
        use libvectors
        use libmap
        implicit none

        integer(c_size_t), intent(in), value :: idx
        real(c_double), intent(in), dimension(4), target :: qvec
        integer(c_size_t), intent(in), value :: nin, nout, nspect, nspin, nlorentz
        complex(c_double_complex), dimension(nlorentz, nspin), intent(out) :: cur
        type(c_ptr), intent(in), value, target :: ff
        ! C++ is row major, while Fortran is column major
        ! This means the moms passed in are transposed
        integer(c_long), dimension(nin), intent(in) :: pids_in
        integer(c_long), dimension(nout), intent(in) :: pids_out
        integer(c_long), dimension(nspect), intent(in) :: pids_spect
        real(c_double), intent(in), dimension(4, nin+nout+nspect) :: moms
        type(fourvector), dimension(nin) :: mom_in
        type(fourvector), dimension(nout) :: mom_out
        type(fourvector), dimension(nspect) :: mom_spect
        type(fourvector) :: qvector
        type(complex_map) :: ff_dict
        integer(c_size_t) :: i
        type(model_holder) :: model

        do i=1,nin
            mom_in(i) = fourvector(moms(1, i), moms(2, i), moms(3, i), moms(4, i))
        enddo

        do i=1,nout
            mom_out(i) = fourvector(moms(1, nin+i), moms(2, nin+i), moms(3, nin+i), moms(4, nin+i))
        enddo

        do i=1,nspect
            mom_spect(i) = fourvector(moms(1, nin+nout+i), moms(2, nin+nout+i), moms(3, nin+nout+i), moms(4, nin+nout+i))
        enddo

        qvector = fourvector(qvec(1), qvec(2), qvec(3), qvec(4))
        ff_dict = complex_map(ff)
        model = models%get_model(idx)
        call model%model_ptr%currents(pids_in, mom_in, nin, pids_out, mom_out, nout, pids_spect, mom_spect, nspect, qvector, ff_dict, cur, nspin, nlorentz)
    end subroutine

    function get_init_wgt(idx, pids_in, pids_spect, pmom, nin, nspect, nproton, nneutron) result(wgt) bind(c, name="GetInitialStateWeight")
        use iso_c_binding
        use libvectors
        implicit none

        integer(c_size_t), intent(in), value :: idx
        integer(c_long), dimension(nin), intent(in) :: pids_in
        integer(c_long), dimension(nspect), intent(in) :: pids_spect
        real(c_double), dimension(4, nin+nspect), intent(in) :: pmom
        type(fourvector), dimension(nin) :: mom_in
        type(fourvector), dimension(nspect) :: mom_spect
        integer(c_size_t), intent(in), value :: nin, nspect, nproton, nneutron
        real(c_double) :: wgt
        integer(c_size_t) :: i
        type(model_holder) :: model

        do i=1,nin
            mom_in(i) = fourvector(pmom(1, i), pmom(2, i), pmom(3, i), pmom(4, i))
        enddo

        do i=1,nspect
            mom_spect(i) = fourvector(pmom(1, nin+i), pmom(2, nin+i), pmom(3, nin+i), pmom(4, nin+i))
        enddo

        model = models%get_model(idx)
        wgt = model%model_ptr%init_wgt(pids_in, mom_in, nin, pids_spect, mom_spect, nspect, nproton, nneutron)
    end function get_init_wgt
end module
