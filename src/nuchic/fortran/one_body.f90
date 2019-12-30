module parameters
    implicit none
    real*8, parameter :: xmp=938.272046d0,xmn=939.56563d0,xm=0.5d0*(xmn+xmp)
    real*8, parameter :: hbarc=197.327053d0,G_F = 1.14e-5*0.1973d0**2
    real*8, parameter :: pi= 4.0d0*atan(1.0d0),c0=1.0d0/16.0d0/pi/pi
    real*8, save :: e,ef,cost,q2,p2,pf2,sint,sina2,wf,wtf
end module parameters



subroutine one_body(xq,w,wt,xk,xp,phi,ee0,theta,ig,sig)
    use parameters
    use dirac_matrices
    implicit none
    integer*4, INTENT IN :: ig!,i
    real*8, INTENT OUT :: sig
    real*8, INTENT IN :: xq, w,wt, xk, xp, ee0, theta, phi
    real*8 :: ek,epf,xmf
    real*8 :: xk_x,xk_y,xk_z
    real*8 :: p_4(4),pp_4(4),q_4(4), cosa
!    ee0      incident electron energy in MeV
!    w        electron energy loss in MeV
!    xq       three-momentum transfer in inverse fm
!    thetae   electron scattering angle
!    xk, xp   initital and final nucleon three-momenta in inverse fm
!    sina2    sin(alpha)**2, alpha being the angle between the 
    e = ee0/hbarc
    wf = w/hbarc
    wtf=wt/hbarc
    ef = e-wf            
    q2 = xq**2
    p2 = xk**2           
    pf2 = xp**2          
    sint = sin(theta)
    cost = sqrt(1.0d0-sint**2)
    if(theta.ge.0.5*pi)cost=-cost
    cosa = 0.d0
    if(xk.ne.0.) cosa=((pf2-p2-q2)/2.0d0/xk/xq)
    sina2 = 1.0d0-cosa**2
    xmf=xm/hbarc
    ek=sqrt(xk**2+xmf**2)
    epf = sqrt(xmf**2 + xp**2)
    !
    !q_4(1)=wtf
    q_4(1)=wtf
    q_4(2:3)=0.0d0
    q_4(4)=xq
    !
    xk_x=xk*sqrt(sina2)*cos(phi)
    xk_y=xk*sqrt(sina2)*sin(phi)      
    xk_z=xk*cosa
    p_4(1)=ek
    p_4(2)=xk_x
    p_4(3)=xk_y
    p_4(4)=xk_z
    pp_4(1)=epf
    pp_4(2)=xk_x
    pp_4(3)=xk_y
    pp_4(4)=xk_z+xq
    call current_init(p_4,pp_4,q_4)
    call define_spinors()
    call sigccc(sig,ig)

    return
end subroutine cc1

subroutine sigccc(sig,ig)
    use xsec
    use dirac_matrices
    implicit none
    integer*4 :: ig
    real*8, parameter :: alpha=1.0d0/137.0d0
    real*8 :: tan2,sig_mott,qm2
    real*8 :: ff1s,ff2s,ff1v,ff2v,ffa,ffp,ges,gms,gev,gmv
    real*8 :: sig,rlp,rtp,rln,rtn,al,at,ff1p,ff1n,ff2p,ff2n

    tan2 = (1.0d0-cost)/(1.0d0+cost)   
    !.....compute sigma_mott [ fm^2 --> mb ]
    if(cost.eq.1.0d0 .or.cost.eq.-1.0d0 )then
        write(6,*)' >>>> warning: divergence of Mott cross section '
        write(6,*)'      theta = ',acos(cost)*180./pi,' deg '
        stop
    else      
        sig_mott = alpha**2/2./(1.0d0-cost)/e**2/tan2   
    end if
    sig_mott = 10.0d0*sig_mott
    qm2 = q2-wf**2
    !qm2 = q2-wtf**2
    al=(-qm2/q2)**2
    at=(qm2/q2/2.0d0+tan2)


    call nform(ig,qm2,ff1s,ff2s,ff1v,ff2v,ffa,ffp,ges,gms,gev,gmv)
    ff1p=0.5d0*(ff1v+ff1s)
    ff2p=0.5d0*(ff2v+ff2s)
    ff1n=0.5d0*(-ff1v+ff1s)
    ff2n=0.5d0*(-ff2v+ff2s)
    call det_Ja(ff1p,ff2p)
    call det_res1b(rlp,rtp)
    call det_Ja(ff1n,ff2n)
    call det_res1b(rln,rtn)
    !write(6,*)'rl,rt',rl,rt
    sig=sig_mott*0.5d0*(al*(rlp+rln)+at*(rtp+rtn))

    return
end subroutine sigccc


subroutine two_body(pp1,ctpp1,p2,ctp2,phip2,p1,ctp1,phip1,ip1,ip2,ie1,ie2,w,qval,r_now)
    use dirac_matrices         
    use mathtool
    implicit none
    real*8, parameter :: lsq=0.71*1.e6,l3=3.5d0*1.e6,xma2=1.1025d0*1.e6
    real*8, parameter :: fstar=2.13d0,eps=20.0d0
    integer*4, INTENT IN :: ie1,ie2,ip1,ip2
    real*8, INTENT IN :: w,ctpp1,p2,ctp2,phip2,p1,ctp1,phip1,qval,pp1
    real*8, INTENT OUT :: r_now
    real*8 :: stpp1,stp1,stp2
    real*8 :: at,bt,vt,par1,par2,pp1,den,jac,arg
    real*8 :: q2,rho,norm,ca5,cv3,gep
    real*8 :: p1_4(4),p2_4(4),pp1_4(4),pp2_4(4),k2_4(4),k1_4(4),q_4(4),pp_4(4)
    real*8 :: k2e_4(4),k1e_4(4)
    real*8 :: r_cc_pi,r_cl_pi,r_ll_pi,r_t_pi,r_tp_pi
    real*8 :: r_cc_del,r_cl_del,r_ll_del,r_t_del,r_tp_del   
    real*8 :: r_cc_int,r_cl_int,r_ll_int,r_t_int,r_tp_int  
    real*8 :: dp1,dp2,delta_w
    real*8 :: tkin_pp1,tkin_pp2, u_pp1,u_pp2
    real*8 :: dir(5),exc(5),r_now(5)


    stpp1=sqrt(1.0d0-ctpp1**2)
    stp1=sqrt(1.0d0-ctp1**2)
    stp2=sqrt(1.0d0-ctp2**2)
    pp1=xpf+ran()*ppmax

    p1_4(1)=sqrt(p1**2+xmn**2)
    p1_4(2)=p1*stp1*cos(phip1)
    p1_4(3)=p1*stp1*sin(phip1)
    p1_4(4)=p1*ctp1
    p2_4(1)=sqrt(p2**2+xmn**2)
    p2_4(2)=p2*stp2*cos(phip2)
    p2_4(3)=p2*stp2*sin(phip2)
    p2_4(4)=p2*ctp2

    !......define constants and ff
    q2=w**2-qval**2
    gep=1.0d0/(1.0d0-q2/lsq)**2 
    !ffgnd= fstar/(1.0d0-q2/lsq)**2/(1.0d0-q2/4.0d0/lsq)*sqrt(3.0d0/2.0d0)
    cv3=fstar/(1.0d0-q2/lsq)**2/(1.0d0-q2/4.0d0/lsq)*sqrt(3.0d0/2.0d0)
    ca5=0.0d0!1.20d0/(1.0d0-q2/xma2)**2/(1.0d0-q2/3.0d0/xma2)*sqrt(3.0d0/2.0d0)
    !   ffgnd=gep/sqrt(1.0d0-q2/(xmn+xmd)**2)/sqrt(1.0d0-q2/l3)
    rho=xpf**3/(1.5d0*pi**2)

    !....Pauli blocking
    if(pp1.lt.xpf) then   
        r_now=0.0d0
        return
    endif        
    !....at this point we can define pp1_4
    pp1_4(1)=sqrt(pp1**2+xmn**2)
    pp1_4(2)=pp1*stpp1
    pp1_4(3)=0.0d0
    pp1_4(4)=pp1*ctpp1


    q_4(2:3)=0.0d0
    q_4(4)=qval
    !...define pp2
    pp2_4(2:4)=p1_4(2:4)+p2_4(2:4)-pp1_4(2:4)+q_4(2:4)
    !....Pauli blocking   
    if(sqrt(sum(pp2_4(2:4)**2)).lt.xpf) then
        r_now=0.0d0
        return
    endif
    !....probably this is not necessary, I need to think about it   
    pp2_4(1)=sqrt(sum(pp2_4(2:4)**2)+xmn**2)
    !...define the argument of the delta-function

    tkin_pp1=pp1_4(1)-xmn
    tkin_pp2=pp2_4(1)-xmn
    u_pp1=0.0d0
    u_pp2=0.0d0
    !if(i_fsi.eq.1.and.(tkin_pp1.lt.kin(npot)).and.(tkin_pp1.gt.kin(1))) call interpolint(kin,pot,npot,tkin_pp1,u_pp1,1)
    !if(i_fsi.eq.1.and.(tkin_pp2.lt.kin(npot)).and.(tkin_pp2.gt.kin(1))) call interpolint(kin,pot,npot,tkin_pp2,u_pp2,1)

    if(i_fg.eq.1) then
        q_4(1)=w-40.0d0
    else
        ! q_4(1)=w-p1_4(1)-p2_4(1)-ep(ie1)+xmn-ep(ie2)+xmn+60.0d0!-u_pp1-u_pp2  
        q_4(1)=w-p1_4(1)-p2_4(1)-ep(ie1)+xmn-ep(ie2)+xmn+120.0d0  
    endif

    !...delta function
    arg=q_4(1)+p1_4(1)+p2_4(1)-pp1_4(1)-pp2_4(1)
    delta_w=fdelta(arg,eps)
    !...define pion momenta
    k1_4(:)=pp1_4(:)-p1_4(:)
    k2_4(:)=q_4(:)-k1_4(:)
    k1e_4(:)=pp2_4(:)-p1_4(:)
    k2e_4(:)=q_4(:)-k1e_4(:)

    !.......currents
    call current_init(p1_4,p2_4,pp1_4,pp2_4,q_4,k1_4,k2_4,1)      
    call define_spinors()
    call det_Jpi(gep)
    call det_JpiJpi(r_cc_pi,r_cl_pi,r_ll_pi,r_t_pi,r_tp_pi)
    call det_JaJb_JcJd(cv3,ca5,np_del,pdel,pot_del)
    call det_JaJc_dir(r_cc_del,r_cl_del,r_ll_del,r_t_del,r_tp_del)
    call det_JpiJaJb(r_cc_int,r_cl_int,r_ll_int,r_t_int,r_tp_int)

    dir(1)=r_cc_pi+2.0d0*(r_cc_del+r_cc_int)
    dir(2)=r_cl_pi+2.0d0*(r_cl_del+r_cl_int)
    dir(3)=r_ll_pi+2.0d0*(r_ll_del+r_ll_int)
    dir(4)=r_t_pi+2.0d0*(r_t_del+r_t_int)
    dir(5)=r_tp_pi+2.0d0*(r_tp_del+r_tp_int)

    call current_init(p1_4,p2_4,pp2_4,pp1_4,q_4,k1e_4,k2e_4,2)
    call det_JaJb_JcJd(cv3,ca5,np_del,pdel,pot_del)
    call det_JaJc_exc(r_cc_del,r_cl_del,r_ll_del,r_t_del,r_tp_del)
    call det_JpiJaJb_exc(r_cc_int,r_cl_int,r_ll_int,r_t_int,r_tp_int)

    exc(1)=2.0d0*(r_cc_del+r_cc_int)
    exc(2)=2.0d0*(r_cl_del+r_cl_int)
    exc(3)=2.0d0*(r_ll_del+r_ll_int)
    exc(4)=2.0d0*(r_t_del+r_t_int)
    exc(5)=2.0d0*(r_tp_del+r_tp_int)


    dp1=PkE(ip1,ie1)
    dp2=PkE(ip2,ie2)


    r_now(:) =dp1*dp2*p1**2*p2**2/(2.0d0*pi)**8*(dir(:)-exc(:))* &
        &      delta_w*pp1**2/rho*dble(xA)/2.0d0/2.0d0 ! /2.0d0 for the electromagnetic piece

    return
end subroutine two_body





