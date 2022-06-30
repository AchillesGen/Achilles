interface

    subroutine get_nucleon_c(nuc, i, part) bind(C, name="GetNucleon")
        use iso_c_binding 
        implicit none

        type(c_ptr) :: nuc
        integer(c_size_t), value :: i
        type(c_ptr) :: part
    end subroutine

    subroutine add_particle_c(nuc, part) bind(C, name="AddParticle")
        use iso_c_binding
        implicit none

        type(c_ptr) :: nuc
        type(c_ptr) :: part
    end subroutine
end interface
