interface

    subroutine event_momentum_c(evt, i, mom) bind(C, name="GetEventMomentum")
        use iso_c_binding
        implicit none

        type(c_ptr) :: evt
        integer(c_size_t), value :: i
        type(c_ptr) :: mom
    end subroutine

    subroutine event_nucleus_c(evt, nuc) bind(C, name="GetEventNucleus")
        use iso_c_binding
        implicit none

        type(c_ptr) :: evt, nuc
    end subroutine
end interface
