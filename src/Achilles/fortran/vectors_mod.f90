module libvectors
    use iso_c_binding

    private
    public :: fourvector, threevector

    include "vectors_cdef.f90"

    type fourvector
        private
        type(c_ptr) :: ptr ! Pointer to the FourVector obj
    contains
        ! Bind some functions to the type for cleaner syntax
        final :: delete_fourvector

        ! Member functions
        procedure :: self => get_ptr4
        procedure :: boost => boost4
        procedure :: boost_vector =>  boost_vector4
        procedure :: dot => dot4
        procedure :: add => add4
        procedure :: sub => sub4
        procedure :: scale => scale4
        procedure :: get => get4
        procedure :: print => print4
        procedure :: to_array => array4
    end type fourvector

    ! This will act as constructor
    interface fourvector
        module procedure create_fourvector
        module procedure copy_constructor4
    end interface

    type threevector
        private
        type(c_ptr) :: ptr ! Pointer to the ThreeVector obj
    contains
        ! Bind some functions to the type for cleaner syntax
        final :: delete_threevector

        ! Function members
        procedure :: self => get_ptr3
        procedure :: dot => dot3
        procedure :: add => add3
        procedure :: sub => sub3
        procedure :: scale => scale3
        procedure :: get => get3
        procedure :: print => print3
        procedure :: copy => new3
    end type threevector

    ! This will act as constructor
    interface threevector
        module procedure create_threevector
        module procedure copy_constructor3
    end interface
    ! end interface

    ! Operator overloads
    interface operator (+)
        module procedure add4, add3
    end interface

    interface operator (-)
        module procedure sub4, sub3
    end interface

    interface operator (*)
        module procedure scale4, scale3
        module procedure scale4_2, scale3_2
        module procedure dot4, dot3
    end interface

    interface operator (/)
        module procedure div4, div3
    end interface

    public :: operator(+), operator(-), operator(*), operator(/)

contains ! Implementation of functions
    ! FourVector
    function create_fourvector(e, x, y, z)
        implicit none
        type(fourvector) :: create_fourvector
        double precision, intent(in) :: e, x, y, z
        create_fourvector%ptr = create_fourvector_c(e, x, y, z)
    end function

    function copy_constructor4(other)
        implicit none
        type(c_ptr), intent(in), target :: other
        type(fourvector) :: copy_constructor4
        copy_constructor4%ptr = other
    end function

    subroutine delete_fourvector(this)
        implicit none
        type(fourvector) :: this
        call delete_fourvector_c(this%ptr)
    end subroutine

    function get_ptr4(this)
        implicit none
        class(fourvector), intent(in) :: this
        type(c_ptr) :: get_ptr4
        get_ptr4 = this%ptr
    end function

    function boost4(this, beta)
        implicit none
        type(fourvector) :: boost4
        class(fourvector), intent(in) :: this
        class(threevector), intent(in) :: beta

        boost4%ptr = boost_c(this%ptr, beta%ptr)
    end function

    function boost_vector4(this)
        implicit none
        type(threevector) :: boost_vector4
        class(fourvector), intent(in) :: this

        boost_vector4%ptr = boost_vector_c(this%ptr)
    end function
    !end function

    function dot4(this, other)
        implicit none
        double precision :: dot4
        class(fourvector), intent(in) :: this, other

        dot4 = dot4_c(this%ptr, other%ptr)
    end function

    function add4(this, other)
        implicit none
        type(fourvector) :: add4
        class(fourvector), intent(in) :: this, other

        add4%ptr = add4_c(this%ptr, other%ptr)
    end function

    function sub4(this, other)
        implicit none
        type(fourvector) :: sub4
        class(fourvector), intent(in) :: this, other

        sub4%ptr = sub4_c(this%ptr, other%ptr)
    end function

    function scale4(this, other)
        implicit none
        type(fourvector) :: scale4
        class(fourvector), intent(in) :: this
        double precision, intent(in) :: other

        scale4%ptr = scale4_c(this%ptr, other)
    end function

    function scale4_2(other, this)
        implicit none
        type(fourvector) :: scale4_2
        class(fourvector), intent(in) :: this
        double precision, intent(in) :: other

        scale4_2%ptr = scale4_c(this%ptr, other)
    end function

    function div4(this, other)
        implicit none
        type(fourvector) :: div4
        class(fourvector), intent(in) :: this
        double precision, intent(in) :: other

        div4%ptr = scale4_c(this%ptr, 1d0/other)
    end function

    double precision function get4(this, idx)
        implicit none
        class(fourvector), intent(in) :: this
        integer :: idx

        get4 = get4_c(this%ptr, idx)
    end function

    subroutine print4(this)
        implicit none
        class(fourvector), intent(in) :: this

        call print4_c(this%ptr)
    end subroutine

    function array4(this)
        implicit none
        class(fourvector), intent(in) :: this
        double precision, dimension(4) :: array4


        array4(1) = get4(this,0)
        array4(2) = get4(this,1)
        array4(3) = get4(this,2)
        array4(4) = get4(this,3)


    end function

    function create_threevector(x, y, z)
        implicit none
        type(threevector) :: create_threevector
        real*8, intent(in) :: x, y, z
        create_threevector%ptr = create_threevector_c(x, y, z)
    end function

    function copy_constructor3(other)
        implicit none
        type(c_ptr), intent(in) :: other
        type(threevector) :: copy_constructor3
        copy_constructor3%ptr = other
    end function

    function new3(this)
        implicit none
        class(threevector), intent(in) :: this
        type(c_ptr) :: new3 
        new3 = new3_c(this%ptr)
    end function

    subroutine delete_threevector(this)
        implicit none
        type(threevector) :: this
        call delete_threevector_c(this%ptr)
    end subroutine

    function get_ptr3(this)
        implicit none
        class(threevector), intent(in) :: this
        type(c_ptr) :: get_ptr3
        get_ptr3 = this%ptr
    end function

    function dot3(this, other)
        implicit none
        double precision :: dot3
        class(threevector), intent(in) :: this, other

        dot3 = dot3_c(this%ptr, other%ptr)
    end function

    function add3(this, other)
        implicit none
        type(threevector) :: add3
        class(threevector), intent(in) :: this, other

        add3%ptr = add3_c(this%ptr, other%ptr)
    end function

    function sub3(this, other)
        implicit none
        type(threevector) :: sub3
        class(threevector), intent(in) :: this, other

        sub3%ptr = sub3_c(this%ptr, other%ptr)
    end function

    function scale3(this, other)
        implicit none
        type(threevector) :: scale3
        class(threevector), intent(in) :: this
        double precision, intent(in) :: other

        scale3%ptr = scale3_c(this%ptr, other)
    end function

    function scale3_2(other, this)
        implicit none
        type(threevector) :: scale3_2
        class(threevector), intent(in) :: this
        double precision, intent(in) :: other

        scale3_2%ptr = scale3_c(this%ptr, other)
    end function

    function div3(this, other)
        implicit none
        type(threevector) :: div3
        class(threevector), intent(in) :: this
        double precision, intent(in) :: other

        div3%ptr = scale3_c(this%ptr, 1d0/other)
    end function

    double precision function get3(this, idx)
        implicit none
        class(threevector), intent(in) :: this
        integer :: idx

        get3 = get3_c(this%ptr, idx)
    end function

    subroutine print3(this)
        implicit none
        class(threevector), intent(in) :: this

        call print3_c(this%ptr)
    end subroutine
end module
