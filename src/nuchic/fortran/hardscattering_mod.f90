module libhardscattering

    private
    public :: onebody, delete_onebody

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

    subroutine delete_onebody() bind(C, name="Delete")
        use quasi_el
        implicit none
        call delete_pke()
    end subroutine

    subroutine onebody_init(this, fname_p, fname_n, fg, iform)
        use quasi_el
        use libutilities
        use dirac_matrices
        implicit none

        class(onebody_calc) :: this
        integer*4 :: fg, nz, na, iform
        real*8 :: kf, mqef, mleptf
        character(len=*), intent(in) :: fname_p, fname_n

        call init(constants)
        call init_pke(fname_p, fname_n, fg, iform)
        mqef=constants%mqe/constants%hbarc
        mleptf=0.0
        call dirac_matrices_in(mqef,mleptf)
    end subroutine

    subroutine onebody_xsec(this, pn_vec, pfn_vec, kl_vec, kfl_vec, e, mom, w, qval, theta, ee, &
                            nZ, nA, kf, results, length)
        use quasi_el
        use libvectors
        implicit none

        class(onebody_calc) :: this
        type(fourvector), intent(in) :: pn_vec, pfn_vec, kl_vec, kfl_vec
        double precision, dimension(4) :: pn, pfn, kl, kfl
        double precision, intent(in) :: w, qval, theta, ee, e, mom, kf
        integer, intent(in) :: nZ, nA
        integer, intent(out) :: length
        double precision, intent(out), dimension(:), pointer :: results

        pn = pn_vec%to_array()
        pfn = pfn_vec%to_array()
        kl = kl_vec%to_array()
        kfl = kfl_vec%to_array()


        length = 2
        allocate(results(length))

        call f_eval(pn, pfn, kl, kfl, e, mom, w, qval, theta, ee, nZ, nA, kf, results(1), results(2))
    end subroutine

    subroutine initialize(name_p, name_n, fg, iform) bind(C, name="InitializeOneBody")
        use iso_c_binding
        use libutilities
        implicit none

        type(c_ptr), intent(in), value :: name_p, name_n
        integer(c_int), intent(in), value :: fg, iform
        character(len=:), allocatable :: fname_p, fname_n

        fname_p = c2fstring(name_p)
        fname_n = c2fstring(name_n)

        call onebody%init(fname_p, fname_n, fg, iform)
    end subroutine

    subroutine cross_section(pn_vec, pfn_vec, kl_vec, kfl_vec, e, mom, w, qval, theta, ee, &
                             nZ, nA, kf, results, length) &
            bind(C, name="CrossSectionOneBody")
        use iso_c_binding
        use libvectors
        implicit none

        type(c_ptr), intent(in) :: pn_vec, pfn_vec
        type(c_ptr), intent(in) :: kl_vec, kfl_vec
        type(fourvector) :: pn, pfn
        type(fourvector) :: kl, kfl
        real(c_double), intent(in), value :: w, qval, theta, ee, e, mom, kf
        integer(c_int), intent(in), value :: nZ, nA
        type(c_ptr), intent(out) :: results
        integer(c_int), intent(out) :: length
        double precision, dimension(:), pointer :: tmp_results

        pn = fourvector(pn_vec)
        pfn = fourvector(pfn_vec)
        kl = fourvector(kl_vec)
        kfl = fourvector(kfl_vec)

        call onebody%xsec(pn, pfn, kl, kfl, e, mom, w, qval, theta, ee, &
                          nZ, nA, kf, tmp_results, length)
        results = c_loc(tmp_results(1))
    end subroutine

    subroutine clean_up(results, length) bind(C, name="CleanUp")
        use iso_c_binding
        implicit none

        type(c_ptr), intent(in) :: results
        integer(c_int), intent(in) :: length
        double precision, dimension(:), pointer :: tmp_results

        call c_f_pointer(results, tmp_results, [length])
        deallocate(tmp_results)
    end subroutine
end module
