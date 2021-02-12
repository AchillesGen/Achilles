module libhardscattering_calculator
    use liblorentz
    implicit none

    private
    complex*16, private, parameter :: czero = (0.0d0,0.0d0)
    complex*16, private, parameter :: ci    = (0.0d0,1.0d0)

    type, abstract, public :: hardscattering_calculator
        private
        contains
            procedure(hs_dirac_matrices), deferred :: dirac_matrices
            procedure(hs_current_init), deferred :: current_init
            procedure(hs_define_spinors), deferred :: define_spinors
            procedure(hs_xsec), deferred :: calculate_xsec
    end type

    abstract interface
        subroutine hs_dirac_matrices(this, xmn_in)
            import
            class(hardscattering_calculator), intent(inout) :: this
            real*8, intent(in) :: xmn_in
        end subroutine

        subroutine hs_current_init(this, mom)
            import
            class(hardscattering_calculator), intent(inout) :: this
            real*8, intent(in) :: mom(:,:)
        end subroutine

        subroutine hs_define_spinors(this)
            import
            class(hardscattering_calculator), intent(inout) :: this
        end subroutine

        subroutine hs_xsec(this, iformfactor, qm2, al, at, sig_mott, sig)
            import
            class(hardscattering_calculator), intent(inout) :: this
            integer, intent(in) :: iformfactor
            real*8, intent(in) :: qm2, al, at, sig_mott
            real*8, intent(out) :: sig(:)
        end subroutine
    end interface
   
    type, extends(hardscattering_calculator), public :: quasielastic
        private
            complex*16 :: up1(2, 4), upp1(2, 4), ubarp1(2, 4), ubarpp1(2, 4)
            real*8 :: p1(4), pp1(4), q(4)
            complex*16 :: J_1(4,4,4)
            real*8 :: xmn
        contains
            procedure :: dirac_matrices => qe_dirac_matrices
            procedure :: current_init => qe_current_init
            procedure :: define_spinors => qe_define_spinors
            procedure :: calculate_xsec => qe_xsec
            procedure :: det_Ja => qe_det_Ja
            procedure :: det_res1b => qe_det_res1b
    end type

    type, extends(hardscattering_calculator), public :: pion_production
        private
            real*8 :: p1(4), pp1(4), q(4), kpi(4), xmn
        contains
            procedure :: dirac_matrices => pion_dirac_matrices
            procedure :: current_init => pion_current_init
            procedure :: define_spinors => pion_define_spinors
            procedure :: calculate_xsec => pion_xsec
            procedure :: det_res1b => pion_det_res1b
    end type

contains

    subroutine qe_dirac_matrices(this, xmn_in)
        use libutilities
        implicit none
        class(quasielastic), intent(inout) :: this
        real*8, intent(in) :: xmn_in
        this%xmn = xmn_in/constants%hbarc
        call lorentz%initialize()
    end subroutine

    subroutine qe_current_init(this, mom)
        use libutilities
        implicit none
        class(quasielastic), intent(inout) :: this
        real*8, intent(in) :: mom(:,:)
        this%p1 = mom(1, :)/constants%hbarc
        this%q = mom(2, :)/constants%hbarc
        this%pp1 = mom(3, :)/constants%hbarc
    end subroutine

    subroutine qe_define_spinors(this)
        implicit none
        class(quasielastic), intent(inout) :: this
        call lorentz%spinor(this%p1, this%xmn, this%up1, this%ubarp1) 
        call lorentz%spinor(this%pp1, this%xmn, this%upp1, this%ubarpp1) 
    end subroutine

    subroutine qe_xsec(this, iformfactor, qm2, al, at, sig_mott, sig)
        implicit none
        class(quasielastic), intent(inout) :: this
        integer, intent(in) :: iformfactor
        real*8, intent(in) :: qm2, al, at, sig_mott
        real*8, intent(out) :: sig(:)
        real*8 :: ff1s, ff2s, ff1v, ff2v, ffa, ffp, ges, gms, gev, gmv
        real*8 :: ff1p, ff2p, ff1n, ff2n
        real*8 :: rlp, rtp, rln, rtn

        ! Load form factors
        call nform(iformfactor, qm2, ff1s, ff2s, ff1v, ff2v, ffa, ffp, ges, gms, gev, gmv)
        ff1p = 0.5d0*(ff1v+ff1s) 
        ff2p = 0.5d0*(ff2v+ff2s) 
        ff1n = 0.5d0*(-ff1v+ff1s) 
        ff2n = 0.5d0*(-ff2v+ff2s) 

        ! Calculate proton xsec
        call this%det_Ja(ff1p, ff2p)
        call this%det_res1b(rlp, rtp)
        sig(1) = sig_mott*0.5d0*(al*rlp+at*rtp)

        ! Calculate proton xsec
        call this%det_Ja(ff1n, ff2n)
        call this%det_res1b(rln, rtn)
        sig(2) = sig_mott*0.5d0*(al*rln+at*rtn)

    end subroutine

    subroutine qe_det_ja(this, f1v, f2v)
        implicit none
        class(quasielastic), intent(inout) :: this
        real*8, intent(in) :: f1v, f2v
        integer*4 :: mu, nu

        do mu=1,4
            this%J_1(:, :, mu) = czero
            do nu=1,4
                this%J_1(:,:,mu)=this%J_1(:,:,mu)+ci*f2v*lorentz%sigma_munu(:,:,mu,nu)&
                    &           *lorentz%g_munu(nu,nu)*this%q(nu)/2.0d0/this%xmn
            enddo
            this%J_1(:,:,mu)=this%J_1(:,:,mu)+f1v*lorentz%gamma_mu(:,:,mu)
        enddo
    end subroutine

    subroutine qe_det_res1b(this, rl, rt)
        implicit none
        class(quasielastic), intent(inout) :: this
        real*8, intent(out) :: rl, rt
        integer*4 :: i1, f1, i
        complex*16 :: J_mu(2,2,4), J_mu_dag(2,2,4)
        real*8 :: res(4, 4)

        res=0.0d0
        do i1=1,2
            do f1=1,2
                do i=1,4
                    J_mu(f1,i1,i)=sum(this%ubarpp1(f1,:)*matmul(this%J_1(:,:,i),this%up1(i1,:))) 
                    J_mu_dag(f1,i1,i)=conjg(J_mu(f1,i1,i))
                    res(i, i) = res(i, i) + J_mu_dag(f1,i1,i)*J_mu(f1,i1,i)
                enddo
            enddo
        enddo
        
        rl=res(1,1)
        rt=res(2,2)+res(3,3)
    endsubroutine

    subroutine pion_dirac_matrices(this, xmn_in)
        use libutilities
        implicit none
        class(pion_production), intent(inout) :: this
        real*8, intent(in) :: xmn_in
        this%xmn = xmn_in/constants%hbarc
        call lorentz%initialize()
    end subroutine

    subroutine pion_current_init(this, mom)
        implicit none
        class(pion_production), intent(inout) :: this
        real*8, intent(in) :: mom(:,:)
        this%p1 = mom(1, :)
        this%q = mom(2, :)
        this%pp1 = mom(3, :)
        this%kpi = mom(4, :)
    end subroutine

    subroutine pion_define_spinors(this)
        implicit none
        class(pion_production), intent(inout) :: this
    end subroutine

    subroutine pion_xsec(this, iformfactor, qm2, al, at, sig_mott, sig)
        implicit none
        class(pion_production), intent(inout) :: this
        integer, intent(in) :: iformfactor
        real*8, intent(in) :: qm2, al, at, sig_mott
        real*8, intent(out) :: sig(:)
        real*8 :: rlp(6), rtp(6), rln(6), rtn(6) 

        call this%det_res1b(rlp, rtp, rln, rtn)

        sig(1:6)=sig_mott*0.5d0*(al*rlp(:)+at*rtp(:))
        sig(7:12)=sig_mott*0.5d0*(al*rln(:)+at*rtn(:))
    end subroutine

    subroutine pion_det_res1b(this, rlp, rtp, rln, rtn)
        use libutilities
        implicit none
        class(pion_production), intent(inout) :: this
        real*8, intent(out) :: rlp(6), rtp(6), rln(6), rtn(6)
        integer*4 :: i1, f1, i, ip
        complex*16 :: J_mu_p(3,6,3,4), J_mu_p_dag(3,6,3,4)
        complex*16 :: J_mu_n(3,6,3,4), J_mu_n_dag(3,6,3,4)
        real*8 :: res_p(4,4,6), res_n(4,4,6)

        ! Outgoing particle channels
        ! n + pi- = 1
        ! n + pi0 = 2
        ! n + pi+ = 3
        ! p + pi- = 4
        ! p + pi0 = 5
        ! p + pi+ = 6

        J_mu_p=(0.0d0,0.0d0)
        J_mu_n=(0.0d0,0.0d0)
        call amplitude(this%q,this%pp1,this%kpi,10,0,1,0,1,J_mu_p(:,3,:,:))
        call amplitude(this%q,this%pp1,this%kpi,10,0,1,1,1,J_mu_p(:,5,:,:))
        call amplitude(this%q,this%pp1,this%kpi,10,0,-1,0,1,J_mu_n(:,2,:,:))
        call amplitude(this%q,this%pp1,this%kpi,10,0,-1,-1,1,J_mu_n(:,4,:,:))
        J_mu_p=J_mu_p/constants%hbarc
        J_mu_n=J_mu_n/constants%hbarc

        res_p = 0d0
        res_n = 0d0
        J_mu_p_dag=0.0d0  
        do i1=1,3
            do ip=1,6
                do f1=1,3
                    do i=1,4
                        J_mu_p_dag(f1,ip,i1,i)=conjg(J_mu_p(f1,ip,i1,i))
                        J_mu_n_dag(f1,ip,i1,i)=conjg(J_mu_n(f1,ip,i1,i))
                        res_p(i,i,ip)=res_p(i,i,ip)+J_mu_p_dag(f1,ip,i1,i)*J_mu_p(f1,ip,i1,i)*this%xmn**2/(this%p1(1)*this%pp1(1))
                        res_n(i,i,ip)=res_n(i,i,ip)+J_mu_n_dag(f1,ip,i1,i)*J_mu_n(f1,ip,i1,i)*this%xmn**2/(this%p1(1)*this%pp1(1))
                    enddo
                enddo
            enddo
        enddo
            
        rlp(:)=res_p(1,1,:)
        rtp(:)=res_p(2,2,:)+res_p(3,3,:)
        rln(:)=res_n(1,1,:)
        rtn(:)=res_n(2,2,:)+res_n(3,3,:)
    end subroutine
        
end module
