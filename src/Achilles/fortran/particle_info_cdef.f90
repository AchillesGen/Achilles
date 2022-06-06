interface
    
    function create_pinfo_c(id) bind(C, name="CreateParticleInfo")
        use iso_c_binding
        implicit none

        type(c_ptr) :: create_pinfo_c
        integer(c_int), value :: id
    end function

    subroutine delete_pinfo_c(self) bind(C, name="DeleteParticleInfo")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: self
    end subroutine

    function name_c(self) bind(C, name="Name")
        use iso_c_binding
        implicit none
        type(c_ptr) :: name_c
        type(c_ptr), value :: self
    end function

    function pid_c(self) bind(C, name="PID")
        use iso_c_binding
        implicit none
        integer(c_int) :: pid_c
        type(c_ptr), value :: self
    end function

    function charge_c(self) bind(C, name="Charge")
        use iso_c_binding
        implicit none
        real(c_double) :: charge_c
        type(c_ptr), value :: self
    end function

    function spin_c(self) bind(C, name="Spin")
        use iso_c_binding
        implicit none
        real(c_double) :: spin_c
        type(c_ptr), value :: self
    end function

    function mass_c(self) bind(C, name="Mass")
        use iso_c_binding
        implicit none
        real(c_double) :: mass_c
        type(c_ptr), value :: self
    end function

    function width_c(self) bind(C, name="Width")
        use iso_c_binding
        implicit none
        real(c_double) :: width_c
        type(c_ptr), value :: self
    end function
end interface
