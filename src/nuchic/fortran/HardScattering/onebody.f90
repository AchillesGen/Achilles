 module quasi_el
    use libutilities
    use libspectral_function
    use libcontract
    use libhardscattering_calculator

    implicit none
    real*8 :: mp,mn,mqe
    real*8 :: hbarc
    real*8, parameter :: pi= 4.0d0*atan(1.0d0)
    
    integer*4, save, private :: fg,iform

    type(spectral_function) :: spectral
    class(hardscattering_calculator), pointer :: qe_calculator
    type(contract) :: contractor

    contains

    function construct_quasielastic()
        class(quasielastic), pointer :: construct_quasielastic
        allocate(construct_quasielastic)
    end function

    subroutine init_pke(fname_pkep,fname_pken,fg_in,iform_in)
        implicit none
        integer*4 :: fg_in,iform_in
        character(len=*) :: fname_pkep,fname_pken
        call init(constants)
        mp = constants%mp
        mn = constants%mn
        mqe = constants%mqe
        hbarc = constants%hbarc
        iform = iform_in
        fg = fg_in
        if(fg.ne.1)then
            spectral = spectral_function(fname_pkep, fname_pken)
        endif
        qe_calculator=>construct_quasielastic()
        call qe_calculator%dirac_matrices(mqe)
    end subroutine

    subroutine f_eval(p_4,pf_4,e,mom,w,qval,thetalept,ee,nZ,nA,kF,fp_o, fn_o)
        implicit none
        real*8, intent(in) :: e, mom
        integer, intent(in) :: nZ, nA 
        real*8, intent(in) :: kF
        real*8 :: p_4(4)
        real*8, intent(in) :: pf_4(4)
        real*8, intent(in) :: w, qval, thetalept, ee
        real*8 :: wt, pkep, pken, xp, xpf, pke
        real*8, intent(out) :: fp_o, fn_o
        real*8 :: sig_p, sig_n, arg, delta_w, norm
        real*8 :: momentum(3, 4), sig(2)

        fp_o=0.0d0      
        fn_o=0.0d0      
        xp=sqrt(sum(p_4(2:4)**2))
        xpf=sqrt(sum(pf_4(2:4)**2))
        p_4(1)=sqrt(xp**2+mqe**2)

        momentum(1, :) = p_4(:)
        momentum(2, 1) = w
        momentum(2, 2:3) = 0
        momentum(2, 4) = qval
        momentum(3, :) = pf_4(:)

        if(fg.ne.1) then
           pkep = spectral%pke(1, e, mom)
           pken = spectral%pke(2, e, mom)
           wt=w-abs(e)+mqe-p_4(1)
           !call cc1(qval/hbarc,w,wt,xp/hbarc,xpf/hbarc,p_4/hbarc,pf_4/hbarc,ee,thetalept,iform,sig_p,sig_n)
           call contractor%xsec(qe_calculator, momentum, wt, ee, thetalept, iform, sig)
           fp_o=pkep*sig(1)/dble(nZ)
           fn_o=pken*sig(2)/dble(nA-nZ)
        else
            if(xp.gt.kF) then
                fp_o=0.0d0   
                fn_o=0.0d0
                return
            endif
            norm=4.0d0*pi/3.0d0*kF**3  
            pke=1.0/norm   
            wt=w-abs(e)
            ! call cc1(qval/hbarc,w,wt,xp/hbarc,xpf/hbarc,p_4/hbarc,pf_4/hbarc,ee,thetalept,iform,sig_p, sig_n)
            call contractor%xsec(qe_calculator, momentum, wt, ee, thetalept, iform, sig)
            fp_o=pke*sig(1)*pf_4(1)/xp/qval
            fn_o=pke*sig(2)*pf_4(1)/xp/qval
        endif
    
        return
    end subroutine f_eval

end module quasi_el
