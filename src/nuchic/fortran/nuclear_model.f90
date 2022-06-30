module nuclear_model
    use iso_c_binding
    use libevent

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
            procedure(nm_states), deferred :: states
            procedure(nm_nspins), deferred :: nspins 
            procedure(nm_fill_nucleus), deferred :: fill_nucleus
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

        subroutine nm_currents(self, event_in, bosons, nbosons, currents, size)
            use libevent
            import model 
            class(model), intent(inout) :: self
            class(event), intent(in) :: event_in
            integer, intent(in) :: nbosons, size
            integer, dimension(nbosons), intent(in) :: bosons
            complex*16, dimension(size), intent(out) :: currents
        end subroutine

        function nm_states(self, info)
            use libprocess_info
            import model 
            class(model), intent(inout) :: self
            type(process_info), intent(inout) :: info
            logical :: nm_states
        end function

        function nm_nspins(self)
            use iso_c_binding
            import model 
            class(model), intent(inout) :: self
            integer(c_size_t) nm_nspins
        end function

        function nm_fill_nucleus(self, evt, xsec, len)
            use libevent
            use iso_c_binding
            import model
            class(model), intent(inout) :: self
            class(event), intent(inout) :: evt
            integer(c_size_t), value :: len
            real(c_double), dimension(len), intent(in) :: xsec
            logical :: nm_fill_nucleus
        end function
    end interface
end module nuclear_model
