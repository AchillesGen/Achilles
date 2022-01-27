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
            procedure :: spinors => onebody_spinors
            procedure :: spectral => onebody_spectral
            procedure :: current => onebody_current
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

    subroutine onebody_spinors(this, qvec, pin, pout)
        use dirac_matrices
        use libvectors
        use libutilities
        use quasi_el
        implicit none
        class(onebody_calc) :: this
        type(fourvector), intent(in) :: qvec, pin, pout
        double precision, dimension(4) :: q, p4, pp4
        double precision :: xp, w, wt, qm2, e

        q = qvec%to_array()
        p4 = pin%to_array()
        pp4 = pout%to_array()
        xp = sqrt(sum(p4(2:4)**2))
        e = p4(1)
        w = q(1)
        p4(1) = sqrt(xp**2+constants%mqe**2)
        !wt = w-abs(e)+constants%mqe-p4(1)
        !q(1) = wt

        call current_init(p4, pp4, q)
        call define_spinors()
    end subroutine

    subroutine onebody_spectral(this, pid, e, xp, spectral)
        use quasi_el
        implicit none

        class(onebody_calc) :: this
        double precision, intent(out) :: spectral
        double precision, intent(in) :: e, xp
        integer, intent(in) :: pid

        if(pid.eq.2212) then
            spectral = pke_p_interp%call(e, xp)
        else
            spectral = pke_n_interp%call(e, xp)
        end if
    end subroutine

    subroutine onebody_current(this, f1, f2, fa, results)
        use dirac_matrices
        use libvectors
        use libutilities
        use quasi_el
        implicit none

        class(onebody_calc) :: this
        complex*16, intent(in) :: f1, f2, fa
        complex*16, intent(inout), dimension(4) :: results

        call det_Ja(f1, f2, fa)
        call det_current(results)
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

    subroutine set_spinors(q_vec, pin_vec, pout_vec) bind(C, name="SetSpinors")
        use iso_c_binding
        use libvectors
        implicit none
        type(c_ptr), intent(in) :: q_vec, pin_vec, pout_vec
        type(fourvector) :: q, pin, pout

        q = fourvector(q_vec)
        pin = fourvector(pin_vec)
        pout = fourvector(pout_vec)
        call onebody%spinors(q, pin, pout)
    end subroutine

    subroutine get_spectral(pid, e, xp, spectral) bind(C, name="GetSpectral")
        use iso_c_binding
        implicit none
        integer(c_int), intent(in), value :: pid
        real(c_double), intent(in), value :: e, xp
        real(c_double), intent(out) :: spectral

        call onebody%spectral(pid, e, xp, spectral)
    end subroutine

    subroutine hadronic_current(f1, f2, fa, results) &
            bind(C, name="HadronicCurrentOneBody")
        use iso_c_binding
        implicit none

        complex(c_double_complex), intent(in), value :: f1, f2, fa
        complex(c_double_complex), intent(inout), dimension(4) :: results

        call onebody%current(f1, f2, fa, results)
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
