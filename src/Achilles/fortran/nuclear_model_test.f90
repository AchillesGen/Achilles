module nuclear_model_test
    use iso_c_binding
    use nuclear_model_interface
    use nuclear_model
    implicit none
    private
    public :: test, register

    type, extends(model) :: test
        contains
            procedure :: init => test_init
!            procedure :: fill_nucleus => test_fill
!            procedure :: nspins => test_spins
!            procedure :: states => test_states
            procedure :: currents => test_currents
            procedure :: ps_name => test_ps
            procedure :: mode => test_mode
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
        call factory%register_model("test", build_test)
    end subroutine

    function test_init(self, filename)
        use libutilities
        class(test), intent(inout) :: self
        character(len=*), intent(in) :: filename
        logical :: test_init
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

    subroutine test_currents(self, pids_in, mom_in, nin, pids_out, mom_out, nout, qvec, ff, len_ff, cur, nspin, nlorentz)
        use iso_c_binding
        use libvectors
        class(test), intent(inout) :: self
        integer(c_size_t), intent(in), value :: nin, nout, len_ff, nspin, nlorentz
        complex(c_double_complex), dimension(len_ff), intent(in) :: ff
        type(fourvector) :: qvec
        type(fourvector), dimension(nin), intent(in) :: mom_in
        type(fourvector), dimension(nout), intent(in) :: mom_out
        integer(c_int), dimension(nin), intent(in) :: pids_in
        integer(c_int), dimension(nout), intent(in) :: pids_out
        complex(c_double_complex), dimension(nspin, nlorentz), intent(out) :: cur

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

!    function test_states(self, info)
!        use libprocess_info
!        use iso_c_binding
!        class(test), intent(inout) :: self
!        type(process_info), intent(inout) :: info
!        logical :: test_states
!        integer(c_size_t) :: in, out
!        integer(c_long), dimension(:), allocatable :: initial, final
!        in = 1
!        out = 1
!        allocate(initial(in))
!        allocate(final(out))
!        initial(1) = 2212
!        final(1) = 2212
!        test_states = .true.
!        call info%add_state(initial, final, in, out)
!        deallocate(initial)
!        deallocate(final)
!    end function

!    function test_spins(self)
!        use iso_c_binding
!        class(test), intent(inout) :: self
!        integer(c_size_t) test_spins
!        test_spins = 4
!    end function

!    function test_fill(self, evt, xsec, len)
!        use libevent
!        use iso_c_binding
!        class(test), intent(inout) :: self
!        class(event), intent(inout) :: evt
!        integer(c_size_t), value :: len
!        real(c_double), dimension(len), intent(in) :: xsec
!        logical :: test_fill
!        test_fill = .true.
!    end function

end module nuclear_model_test
