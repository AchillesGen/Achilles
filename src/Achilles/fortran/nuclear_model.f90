module nuclear_model
    use iso_c_binding
    use libevent
    implicit none

    private
    public :: model

    type, abstract :: model 
        private 
        contains
            ! Member functions
            procedure(nm_init), deferred :: init
            procedure(nm_cleanup), deferred :: cleanup
            procedure(nm_mode), deferred :: mode
            procedure :: frame => nm_frame
            procedure(nm_name), deferred, nopass :: model_name
            procedure(nm_psname), deferred :: ps_name
            procedure(nm_currents), deferred :: currents 
            procedure(nm_init_wgt), deferred :: init_wgt
            procedure(nm_inspirehep), deferred, nopass :: inspirehep
    end type model

    abstract interface
        function nm_init(self, filename, params)
            use libmap
            import model
            class(model), intent(inout) :: self
            character(len=*), intent(in) :: filename
            type(map), intent(in) :: params
            logical :: nm_init
        end function

        subroutine nm_cleanup(self)
            import model
            class(model), intent(inout) :: self
        end subroutine

        function nm_mode(self)
            import model 
            class(model), intent(inout) :: self
            integer :: nm_mode 
        end function

        function nm_name()
            character(len=:), allocatable :: nm_name
        end function

        function nm_inspirehep()
            character(len=:), allocatable :: nm_inspirehep
        end function

        function nm_psname(self)
            import model 
            class(model), intent(inout) :: self
            character(len=:), allocatable :: nm_psname
        end function

        subroutine nm_currents(self, pids_in, mom_in, nin, pids_out, mom_out, nout, pids_spect, mom_spect, nspect, qvec, ff, cur, nspin, nlorentz)
            use libvectors
            use libmap
            use iso_c_binding
            import model 
            class(model), intent(inout) :: self
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
        end subroutine

        function nm_init_wgt(self, pids_in, mom_in, nin, pids_spect, mom_spect, nspect, nproton, nneutron) result(wgt)
            use libvectors
            use iso_c_binding
            import model 
            class(model), intent(inout) :: self
            integer(c_size_t), intent(in), value :: nin, nspect, nproton, nneutron
            integer(c_long), dimension(nin), intent(in) :: pids_in
            integer(c_long), dimension(nspect), intent(in) :: pids_spect
            type(fourvector), dimension(nin), intent(in) :: mom_in
            type(fourvector), dimension(nspect), intent(in) :: mom_spect
            real(c_double) :: wgt
        end function

    end interface
contains
    function nm_frame(self)
        class(model), intent(inout) :: self
        integer :: nm_frame
        nm_frame = 0
    end function
end module nuclear_model
