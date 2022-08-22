interface

    function process_model_c(info) bind(C, name="ProcessModel")
        use iso_c_binding 
        implicit none

        type(c_ptr) :: info, process_model_c
    end function

    subroutine process_ids_c(info, ids, len) bind(C, name="ProcessIDs")
        use iso_c_binding
        implicit none

        type(c_ptr) :: info
        integer(c_long), dimension(:), intent(out) :: ids
        integer(c_size_t), intent(out) :: len
    end subroutine

    function process_multiplicity_c(info) bind(C, name="ProcessMultiplicity")
        use iso_c_binding
        implicit none

        type(c_ptr) :: info
        integer(c_size_t) :: process_multiplicity_c
    end function

    subroutine process_masses_c(info, masses, len) bind(C, name="ProcessMasses")
        use iso_c_binding
        implicit none

        type(c_ptr) :: info
        integer(c_size_t), value :: len
        real(c_double), dimension(:) :: masses
    end subroutine

    subroutine process_add_state_c(info, initial, final, in , out) bind(C, name="ProcessAddState")
        use iso_c_binding
        implicit none

        type(c_ptr), intent(inout) :: info
        integer(c_size_t), intent(in), value :: in, out
        integer(c_long), dimension(:), intent(in) :: initial(in), final(out)
    end subroutine
end interface
