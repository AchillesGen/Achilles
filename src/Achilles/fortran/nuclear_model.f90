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
            procedure(nm_psname), deferred :: ps_name
            procedure(nm_currents), deferred :: currents 
            procedure(nm_init_wgt), deferred :: init_wgt
    end type model

    abstract interface
        function nm_init(self, filename)
            import model
            class(model), intent(inout) :: self
            character(len=*), intent(in) :: filename
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

        function nm_psname(self)
            import model 
            class(model), intent(inout) :: self
            character(len=:), allocatable :: nm_psname
        end function

        subroutine nm_currents(self, pids_in, mom_in, nin, pids_out, mom_out, nout, qvec, ff, cur, nspin, nlorentz)
            use libvectors
            use libmap
            use iso_c_binding
            import model 
            class(model), intent(inout) :: self
            integer(c_size_t), intent(in), value :: nin, nout, nspin, nlorentz
            type(complex_map), intent(in) :: ff
            type(fourvector) :: qvec
            integer(c_long), dimension(nin), intent(in) :: pids_in
            integer(c_long), dimension(nout), intent(in) :: pids_out
            type(fourvector), dimension(nin), intent(in) :: mom_in
            type(fourvector), dimension(nout), intent(in) :: mom_out
            complex(c_double_complex), dimension(nlorentz, nspin), intent(out) :: cur
        end subroutine

        function nm_init_wgt(self, pids, moms, nin, nproton, nneutron) result(wgt)
            use libvectors
            use iso_c_binding
            import model 
            class(model), intent(inout) :: self
            integer(c_size_t), intent(in), value :: nin, nproton, nneutron
            integer(c_long), dimension(nin), intent(in) :: pids
            type(fourvector), dimension(nin), intent(in) :: moms
            real(c_double) :: wgt
        end function

    end interface
end module nuclear_model
