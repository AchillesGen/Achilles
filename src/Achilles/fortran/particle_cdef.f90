interface
    function create_particle_c(id, momentum, position, istatus) &
        bind(C, name="CreateParticle")
        use iso_c_binding
        implicit none

        type(c_ptr) :: create_particle_c
        integer(c_int), intent(in), value :: id
        type(c_ptr), intent(in), value :: momentum
        type(c_ptr), intent(in), value :: position
        integer(c_int), intent(in), value :: istatus
    end function

    function copy_particle_c(other) bind(C, name="CopyParticle")
        use iso_c_binding
        implicit none
        type(c_ptr), intent(in), value :: other
        type(c_ptr) :: copy_particle_c
    end function

    subroutine delete_particle_c(self) bind(C, name="DeleteParticle")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: self
    end subroutine

    function particle_status_c(self) bind(C, name="GetParticleStatus")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: self
        integer(c_int) :: particle_status_c
    end function

    function particle_info_c(self) bind(C, name="GetParticleInfo")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: self
        type(c_ptr) :: particle_info_c
    end function

    function particle_momentum_c(self) bind(C, name="GetParticleMomentum")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: self
        type(c_ptr) :: particle_momentum_c
    end function

    function particle_position_c(self) bind(C, name="GetParticlePosition")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: self
        type(c_ptr) :: particle_position_c
    end function

    function set_particle_momentum_c(self, momentum) bind(C, name="SetParticleMomentum")
        use iso_c_binding
        implicit none
        type(c_ptr), intent(in), value :: self
        type(c_ptr), intent(in), value :: momentum
        type(c_ptr) :: set_particle_momentum_c
    end function

    function set_particle_position_c(self, position) bind(C, name="SetParticlePosition")
        use iso_c_binding
        implicit none
        type(c_ptr), intent(in), value :: self
        type(c_ptr), intent(in), value :: position
        type(c_ptr) :: set_particle_position_c
    end function
end interface
