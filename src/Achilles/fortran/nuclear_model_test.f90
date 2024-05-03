module nuclear_model_test
    use iso_c_binding
    use nuclear_model_interface
    use nuclear_model
    use libspectral_function
    implicit none
    private
    public :: test, register

    type(spectral_function) :: spectral_p, spectral_n

    type, extends(model) :: test
        contains
            procedure :: init => test_init
            procedure :: currents => test_currents
            procedure, nopass :: model_name => test_name
            procedure :: ps_name => test_ps
            procedure :: mode => test_mode
            procedure :: init_wgt => test_init_wgt
            procedure :: cleanup => test_cleanup
    end type

contains

    subroutine expected_version(version) bind(C, name="ExpectedVersion")
        integer(c_int), dimension(3), intent(inout) :: version
        version(1) = 1
        version(2) = 0
        version(3) = 0
    end subroutine

    subroutine register() bind(C, name="Register")
        type(test) :: model
        call factory%register_model(model%model_name(), build_test)
    end subroutine

    function test_init(self, filename, params)
        use libutilities
        class(test), intent(inout) :: self
        integer :: ios, i
        character(len=*), intent(in) :: filename
        type(map), intent(in) :: params
        character(len=200) :: string
        integer, parameter :: read_unit = 99
        logical :: test_init

        open(unit=read_unit, file=trim(filename), iostat=ios)
        if( ios /= 0 ) then
            test_init = .false.
            return
        endif

        read(read_unit, '(A)', iostat=ios) string
        spectral_p = spectral_function(string)

        read(read_unit, '(A)', iostat=ios) string
        spectral_n = spectral_function(string)

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
        test_mode = 1
    end function

    function test_name()
        character(len=:), allocatable :: test_name
        test_name = "test"
    end function

    function test_ps(self)
        class(test), intent(inout) :: self
        character(len=:), allocatable :: test_ps
        test_ps = "QESpectral"
    end function

    subroutine test_currents(self, pids_in, mom_in, nin, pids_out, mom_out, nout, pids_spect, mom_spect, nspect, qvec, ff, cur, nspin, nlorentz)
        use iso_c_binding
        use libvectors
        use libmap
        class(test), intent(inout) :: self
        integer(c_size_t), intent(in), value :: nin, nout, nspect, nspin, nlorentz
        type(complex_map), intent(in) :: ff
        type(fourvector) :: qvec
        integer(c_long), dimension(nin), intent(in) :: pids_in
        integer(c_long), dimension(nout), intent(in) :: pids_out
        integer(c_long), dimension(nspect), intent(in) :: pids_spect
        type(fourvector), dimension(nin), intent(in) :: mom_in
        type(fourvector), dimension(nout), intent(in) :: mom_out
        type(fourvector), dimension(nspect), intent(in) :: mom_spect
        complex(c_double_complex), dimension(nlorentz, nspin), intent(out) :: cur

        cur(1, 1) = 1
        cur(1, 2) = 2
        cur(1, 3) = 3
        cur(1, 4) = 4
        cur(2, 1) = 5
        cur(2, 2) = 6
        cur(2, 3) = 7
        cur(2, 4) = 8
        cur(3, 1) = 9
        cur(3, 2) = 10
        cur(3, 3) = 11
        cur(3, 4) = 12
        cur(4, 1) = 13
        cur(4, 2) = 14
        cur(4, 3) = 15
        cur(4, 4) = 16
    end subroutine

    function test_init_wgt(self, pids_in, mom_in, nin, pids_spect, mom_spect, nspect, nproton, nneutron) result(wgt)
        use iso_c_binding
        use libvectors
        use libutilities

        class(test), intent(inout) :: self
        integer(c_long), dimension(nin), intent(in) :: pids_in
        type(fourvector), dimension(nin), intent(in) :: mom_in
        type(fourvector), dimension(nspect), intent(in) :: mom_spect
        integer(c_long), dimension(nspect), intent(in) :: pids_spect
        integer(c_size_t), intent(in), value :: nin, nspect, nproton, nneutron
        real(c_double) :: wgt

        wgt = 1
    end function test_init_wgt
end module nuclear_model_test
