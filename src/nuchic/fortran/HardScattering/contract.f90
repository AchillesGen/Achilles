module libcontract
    implicit none

    private
    real*8, parameter :: alpha=1.0d0/137.0d0

    type, public :: contract
        private
            real*8 :: el_energy, el_energyf, q2, cost, omegaf, omegatf
        contains
            procedure :: xsec => contract_xsec
            procedure :: evaluate_xsec => contract_evaluate_xsec
    end type contract

contains

    subroutine contract_xsec(this, mode, mom, omegat, el_energy0, theta, iformfactor, sig) 
        use libutilities
        use libhardscattering_calculator
        implicit none
        class(contract) :: this
        class(hardscattering_calculator), pointer :: mode
        real*8, intent(in) :: mom(:, :), omegat, el_energy0, theta 
        integer*4, intent(in) :: iformfactor
        real*8, intent(out) :: sig(:)

        ! momentum labels: p_in, q, p1_out, ..., pn_out
        this%el_energy = el_energy0/constants%hbarc
        this%omegaf = mom(2, 1)/constants%hbarc
        this%omegatf = omegat/constants%hbarc
        this%el_energyf = this%el_energy - this%omegaf
        this%q2 = sum(mom(2, 2:4)**2)/constants%hbarc/constants%hbarc
        this%cost = cos(theta)

        call mode%current_init(mom)
        call mode%define_spinors()
        call this%evaluate_xsec(mode, iformfactor, sig)
    end subroutine contract_xsec

    subroutine contract_evaluate_xsec(this, mode, iformfactor, sig)
        use libhardscattering_calculator
        implicit none
        class(contract) :: this
        class(hardscattering_calculator), pointer :: mode
        integer*4, intent(in) :: iformfactor
        real*8, intent(out) :: sig(:)
        real*8 :: sig_mott, tan2, qm2, al, at

        tan2 = (1.0d0 - this%cost)/(1.0d0 + this%cost)
        sig_mott = 10.0d0*alpha**2/2./(1.0d0-this%cost)/this%el_energy**2/tan2
        qm2 = this%q2 - this%omegaf**2
        al = (-qm2/this%q2)**2
        at = (qm2/this%q2/2.0d0 + tan2)

        call mode%calculate_xsec(iformfactor, qm2, al, at, sig_mott, sig)
    end subroutine

end module
