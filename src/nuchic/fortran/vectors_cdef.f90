interface

    function create_fourvector_c(x, y, z, e) bind(C, name="CreateFourVector")
        use iso_c_binding
        implicit none

        type(c_ptr) :: create_fourvector_c
        real(c_double), value :: x
        real(c_double), value :: y 
        real(c_double), value :: z
        real(c_double), value :: e 
    end function

    subroutine delete_fourvector_c(self) bind(C, name="DeleteFourVector")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: self
    end subroutine

    function boost_c(self, beta) bind(C, name="Boost")
        use iso_c_binding
        implicit none
        type(c_ptr) :: boost_c
        type(c_ptr), intent(in), value :: self
        type(c_ptr), intent(in), value :: beta
    end function

    function boost_vector_c(self) bind(C, name="BoostVector")
        use iso_c_binding
        implicit none
        type(c_ptr) :: boost_vector_c
        type(c_ptr), intent(in), value :: self
    end function

    function dot4_c(self, other) bind(C, name="Dot4")
        use iso_c_binding
        implicit none
        real(c_double) :: dot4_c
        type(c_ptr), intent(in), value :: self
        type(c_ptr), intent(in), value :: other
    end function

    function add4_c(self, other) bind(C, name="Add4")
        use iso_c_binding
        implicit none
        type(c_ptr) :: add4_c
        type(c_ptr), intent(in), value :: self
        type(c_ptr), intent(in), value :: other
    end function

    function sub4_c(self, other) bind(C, name="Sub4")
        use iso_c_binding
        implicit none
        type(c_ptr) :: sub4_c
        type(c_ptr), intent(in), value :: self
        type(c_ptr), intent(in), value :: other
    end function

    function scale4_c(self, other) bind(C, name="Scale4")
        use iso_c_binding
        implicit none
        type(c_ptr) :: scale4_c
        type(c_ptr), intent(in), value :: self
        real(c_double), intent(in), value :: other
    end function

    function get4_c(self, other) bind(C, name="Get4")
        use iso_c_binding
        implicit none
        real(c_double) :: get4_c
        type(c_ptr), intent(in), value :: self
        integer(c_int), intent(in), value :: other
    end function

    subroutine print4_c(self) bind(C, name="Print4")
        use iso_c_binding
        implicit none
        type(c_ptr), intent(in), value :: self
    end subroutine

    function create_threevector_c(x, y, z) bind(C, name="CreateThreeVector")
        use iso_c_binding
        implicit none

        type(c_ptr) :: create_threevector_c
        real(c_double), value :: x, y, z 
    end function

    function new3_c(this) bind(C, name="New3")
        use iso_c_binding
        implicit none
        type(c_ptr) :: new3_c
        type(c_ptr), intent(in), value :: this
    end function

    subroutine delete_threevector_c(self) bind(C, name="DeleteThreeVector")
        use iso_c_binding
        implicit none
        type(c_ptr), value :: self
    end subroutine

    function dot3_c(self, other) bind(C, name="Dot3")
        use iso_c_binding
        implicit none
        real(c_double) :: dot3_c
        type(c_ptr), intent(in), value :: self
        type(c_ptr), intent(in), value :: other
    end function

    function add3_c(self, other) bind(C, name="Add3")
        use iso_c_binding
        implicit none
        type(c_ptr) :: add3_c
        type(c_ptr), intent(in), value :: self
        type(c_ptr), intent(in), value :: other
    end function

    function sub3_c(self, other) bind(C, name="Sub3")
        use iso_c_binding
        implicit none
        type(c_ptr) :: sub3_c
        type(c_ptr), intent(in), value :: self
        type(c_ptr), intent(in), value :: other
    end function

    function scale3_c(self, other) bind(C, name="Scale3")
        use iso_c_binding
        implicit none
        type(c_ptr) :: scale3_c
        type(c_ptr), intent(in), value :: self
        real(c_double), intent(in), value :: other
    end function

    function get3_c(self, other) bind(C, name="Get3")
        use iso_c_binding
        implicit none
        real(c_double) :: get3_c
        type(c_ptr), intent(in), value :: self
        integer(c_int), intent(in), value :: other
    end function

    subroutine print3_c(self) bind(C, name="Print3")
        use iso_c_binding
        implicit none
        type(c_ptr), intent(in), value :: self
    end subroutine

end interface
