module libhardscattering

    private
    public :: onebody

!    type, abstract hardscattering
!        contains
!            procedure(gen_xsec), deferred :: cross_section
!            procedure(gen_initial), deferred :: initial
!    end type

    type onebody_calc
        contains
            procedure :: init => onebody_init
            procedure :: xsec => onebody_xsec
    end type

    type(onebody_calc) :: onebody

contains

    subroutine onebody_init(this, fname_p, fname_n, fg, nz, na, kf, iform)
        use quasi_el
        use libutilities
        implicit none

        class(onebody_calc) :: this
        integer*4 :: fg, nz, na, iform
        real*8 :: kf
        character(len=*), intent(in) :: fname_p, fname_n

        call init(constants)
        call init_pke(fname_p, fname_n, fg, nz, na, kf, iform)
    end subroutine

    function onebody_xsec(this, inucleon, pn_vec, pfn_vec, e, mom, w, qval, theta, ee)
        use quasi_el
        use libvectors
        implicit none

        class(onebody_calc) :: this
        integer*4 :: inucleon
        type(fourvector), intent(in) :: pn_vec, pfn_vec
        double precision, dimension(4) :: pn, pfn
        double precision :: w, qval, theta, ee, e, mom
        double precision :: onebody_xsec

        pn = pn_vec%to_array()
        pfn = pfn_vec%to_array()

        call f_eval(inucleon, pn, pfn, e, mom, w, qval, theta, ee, onebody_xsec)
    end function

    subroutine initialize(name_p, name_n, fg, nz, na, kf, iform) &
            bind(C, name="InitializeOneBody")
        use iso_c_binding
        use libutilities
        implicit none

        type(c_ptr), intent(in), value :: name_p, name_n
        integer(c_int), intent(in), value :: fg, nz, na, iform
        real(c_double), intent(in), value :: kf
        character(len=:), allocatable :: fname_p, fname_n

        fname_p = c2fstring(name_p)
        fname_n = c2fstring(name_n)

        call onebody%init(fname_p, fname_n, fg, nz, na, kf, iform)
    end subroutine

    function cross_section(inucleon, pn_vec, pfn_vec, e, mom, w, qval, theta, ee) &
            bind(C, name="CrossSectionOneBody")
        use iso_c_binding
        use libvectors
        implicit none

        integer(c_int), intent(in), value :: inucleon
        type(c_ptr), intent(in) :: pn_vec, pfn_vec
        type(fourvector) :: pn, pfn
        real(c_double), intent(in), value :: w, qval, theta, ee, e, mom
        real(c_double) :: cross_section

        pn = fourvector(pn_vec)
        pfn = fourvector(pfn_vec)

        cross_section = onebody%xsec(inucleon, pn, pfn, e, mom, w, qval, theta, ee)
    end function

end module
