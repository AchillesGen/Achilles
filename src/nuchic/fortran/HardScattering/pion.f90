module libpion_production
    use libutilities
    use libspectral_function
    use libcontract
    use libhardscattering_calculator

    implicit none
    real*8 :: mp,mn,mqe
    real*8 :: mpip,mpi0,mpi
    real*8 :: hbarc
    real*8, parameter :: pi = 4.0*atan(1.0d0)

    integer*4, save, private :: fg,iform

    type(spectral_function) :: spectral
    class(hardscattering_calculator), pointer :: pion_calculator
    type(contract) :: contractor

    contains

    function construct_pion()
        class(pion_production), pointer :: construct_pion
        allocate(construct_pion)
    end function

    subroutine init_pke(fname_pkep, fname_pken, fg_in, iform_in)
        use libpartinfo
        implicit none
        integer*4 :: fg_in,i,j,iform_in
        real*8 :: he_p,he_n
        real*8 :: norm,pmax, hp
        real*8 :: dpe, ee, pp, pke
        character(len=*) :: fname_pkep,fname_pken
        type(pinfo) :: pip, pi0, proton, neutron

        call init(constants)
        pip = pinfo(211)
        pi0 = pinfo(111)
        proton = pinfo(2212)
        neutron = pinfo(2112)
        mpip = pip%mass()
        mpi0 = pi0%mass()
        mpi = (2.0*mpip+mpi0)/3.0
        mp = proton%mass()
        mn = neutron%mass()
        mqe = (mp+mn)/2.0
        hbarc = constants%hbarc
        iform= iform_in
        fg = fg_in
        if(fg.ne.1)then
            spectral = spectral_function(fname_pkep, fname_pken)
        endif
        pion_calculator=>construct_pion()
        call pion_calculator%dirac_matrices(mqe)
    end subroutine

    subroutine f_eval(p_4,pf_4,k_4,e,w,qval,thetalept,ee,nZ,nA,kf,fp_o,fn_o)
        implicit none

        ! Move these parameters to the phase space generation
        !real*8, parameter :: k_max=10.0e3, e_min=1076.957d0
        !
        integer*4, intent(in) :: nA,nZ
        real*8, intent(in) :: e,w,kf,thetalept,ee,qval
        real*8 :: p_4(4), pf_4(4), k_4(4)
        real*8, intent(out) :: fp_o(6), fn_o(6)
        real*8 :: wt, xp, xpf, norm, pke, pkep, pken
        real*8 :: sig(12), momentum(4, 4)

        fp_o = 0.d0
        fn_o = 0.d0
        xp = sqrt(sum(p_4(2:4)**2))
        xpf = sqrt(sum(pf_4(2:4)**2))
        p_4(1) = sqrt(xp**2+mqe**2)

        momentum(1, :) = p_4(:)
        momentum(2, 1) = w
        momentum(2, 2:3) = 0
        momentum(2, 4) = qval
        momentum(3, :) = pf_4(:)
        momentum(4, :) = k_4(:)

        if(fg.ne.1) then
            pkep = spectral%pke(1, e, xp)
            pken = spectral%pke(2, e, xp)
            wt = w - abs(e) + mqe - p_4(1)
            ! call cc1(qval/hbarc,w,p_4,pf_4,k_4,ee,thetalept,iform,sig)
            call contractor%xsec(pion_calculator, momentum, wt, ee, thetalept, iform, sig)
            fp_o(:) = pkep*sig(1:6)/dble(nZ)
            fn_o(:) = pken*sig(7:12)/dble(nA-nZ)
        else
            if(xp.gt.kf) then
                fp_o = 0.0d0
                fn_o = 0.0d0
                return
            endif
            norm = 4.0d0*pi/3.0d0*kf**3
            pke=1.0/norm
            wt=w-abs(e)
            ! call cc1(qval/hbarc,w,p_4,pf_4,k_4,ee,thetalept,iform,sig)
            call contractor%xsec(pion_calculator, momentum, wt, ee, thetalept, iform, sig)
            fp_o(:) = pke*sig(1:6)*pf_4(1)/xp/qval
            fn_o(:) = pke*sig(7:12)*pf_4(1)/xp/qval
        endif

        return
    end subroutine

end module
