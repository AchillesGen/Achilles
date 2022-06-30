interface

    subroutine event_momentum_c(evt, i, mom) bind(C, name="GetEventMomentum")
        use iso_c_binding
        implicit none

        type(c_ptr) :: evt
        integer(c_size_t), value :: i
        type(c_ptr) :: mom
    end subroutine

    subroutine event_matrix_wgt_c(evt, i, wgt) bind(C, name="SetEventMatrixWgt")
        use iso_c_binding
        implicit none

        type(c_ptr) :: evt
        integer(c_size_t), value :: i
        real(c_double), intent(in), value :: wgt
    end subroutine

    function event_total_cross_section_c(evt) bind(C, name="EventTotalCrossSection")
        use iso_c_binding
        implicit none

        type(c_ptr) :: evt
        real(c_double) :: event_total_cross_section_c
    end function

    function event_select_nucleon_c(evt) bind(C, name="EventSelectNucleon")
        use iso_c_binding
        implicit none

        type(c_ptr) :: evt
        integer(c_size_t) :: event_select_nucleon_c
    end function

    subroutine event_nucleus_c(evt, nuc) bind(C, name="GetEventNucleus")
        use iso_c_binding
        implicit none

        type(c_ptr) :: evt, nuc
    end subroutine

    subroutine set_event_mewgt_c(evt, wgt) bind(C, name="SetEventMEWgt")
        use iso_c_binding
        implicit none

        type(c_ptr) :: evt
        real(c_double) :: wgt
    end subroutine
end interface
