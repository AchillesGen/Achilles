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
            procedure :: tensor => onebody_tensor
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
        real*8 :: kf, mqef
        character(len=*), intent(in) :: fname_p, fname_n

        call init(constants)
        call init_pke(fname_p, fname_n, fg, iform)
        mqef=constants%mqe/constants%hbarc
        call dirac_matrices_in(mqef)
    end subroutine

    subroutine onebody_xsec(this, pn_vec, pfn_vec, e, mom, w, qval, theta, ee, &
                            nZ, nA, kf, results, length)
        use quasi_el
        use libvectors
        implicit none

        class(onebody_calc) :: this
        type(fourvector), intent(in) :: pn_vec, pfn_vec
        double precision, dimension(4) :: pn, pfn
        double precision, intent(in) :: w, qval, theta, ee, e, mom, kf
        integer, intent(in) :: nZ, nA
        integer, intent(out) :: length
        double precision, intent(out), dimension(:), pointer :: results

        pn = pn_vec%to_array()
        pfn = pfn_vec%to_array()

        length = 2
        allocate(results(length))

        call f_eval(pn, pfn, e, mom, w, qval, theta, ee, nZ, nA, kf, results(1), results(2))
    end subroutine

    subroutine onebody_tensor(this, qvec, pin, pout, results)
        use dirac_matrices
        use libvectors
        use libutilities
        use quasi_el
        implicit none

        class(onebody_calc) :: this
        type(fourvector), intent(in) :: qvec, pin, pout
        double precision, dimension(4) :: q, p4, pp4
        double precision :: f1v, f2v
        complex*16, intent(inout), dimension(16) :: results
        double precision :: ff1s, ff2s, ff1v, ff2v, ffa, ffp, ges, gms, gev, gmv
        double precision :: xp, wt, e, w, qm2
        double precision :: sf

        q = qvec%to_array()
        p4 = pin%to_array()
        pp4 = pout%to_array()
        xp = sqrt(sum(p4(2:4)**2))
        e = p4(1)
        sf = pke_p_interp%call(e, xp) 
        ! print*, 'Spectral: ', e, xp, sf
        w = q(1)
        p4(1) = sqrt(xp**2+constants%mqe**2)
        wt = w-abs(e)+constants%mqe-p4(1)
        q(1) = wt

        q = q/constants%hbarc
        p4 = p4/constants%hbarc
        pp4 = pp4/constants%hbarc

        qm2 = sum(q(2:4)**2) - (w/constants%hbarc)**2

        call current_init(p4, pp4, q)
        call define_spinors()
        call nform(2, qm2, ff1s, ff2s, ff1v, ff2v, ffa, ffp, ges, gms, gev, gmv)
        f1v = 0.5*(ff1v+ff1s)
        f2v = 0.5*(ff2v+ff2s)

        call det_Ja(f1v, f2v, 0d0)
        call det_current(results)
        results = results*sf/6
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

    subroutine cross_section(pn_vec, pfn_vec, e, mom, w, qval, theta, ee, &
                             nZ, nA, kf, results, length) &
            bind(C, name="CrossSectionOneBody")
        use iso_c_binding
        use libvectors
        implicit none

        type(c_ptr), intent(in) :: pn_vec, pfn_vec
        type(fourvector) :: pn, pfn
        real(c_double), intent(in), value :: w, qval, theta, ee, e, mom, kf
        integer(c_int), intent(in), value :: nZ, nA
        type(c_ptr), intent(out) :: results
        integer(c_int), intent(out) :: length
        double precision, dimension(:), pointer :: tmp_results

        pn = fourvector(pn_vec)
        pfn = fourvector(pfn_vec)

        call onebody%xsec(pn, pfn, e, mom, w, qval, theta, ee, &
                          nZ, nA, kf, tmp_results, length)
        results = c_loc(tmp_results(1))
    end subroutine

    subroutine hadronic_current(q_vec, pin_vec, pout_vec, results) &
            bind(C, name="HadronicTensorOneBody")
        use iso_c_binding
        use libvectors
        implicit none

        type(c_ptr), intent(in) :: q_vec, pin_vec, pout_vec
        type(fourvector) :: q, pin, pout
        complex(c_double_complex), intent(inout), dimension(16) :: results
        integer(c_int) :: length

        q = fourvector(q_vec)
        pin = fourvector(pin_vec)
        pout = fourvector(pout_vec)
        call onebody%tensor(q, pin, pout, results)
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
