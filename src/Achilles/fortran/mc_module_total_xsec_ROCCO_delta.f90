module mc_module
   implicit none
   integer*4, private, save :: nev,xA,i_fg,i_fsi,np,ne,nwlk,npot,np_del
   integer*4, private, parameter :: neq=10000,nvoid=10
   real*8, private, save ::  xpf
   real*8, private, save:: xmpi,xmd,xmn,norm,thetalept
   real*8, private, parameter :: pi=acos(-1.0d0),hbarc=197.327053d0,xmp=938.0d0,ppmax=1.0d0*1.e3
   real*8, private,parameter :: alpha=1.0d0/137.0d0
   real*8, private, allocatable :: pv(:),dp(:),ep(:),Pke(:,:)
   real*8, private, allocatable :: kin(:),pot(:),pdel(:),pot_del(:)
   integer*8, private, allocatable, save :: irn(:)
contains

subroutine mc_init(i_fg_in,pwia,i_fsi_in,irn_in,nev_in,nwlk_in,xpf_in,thetalept_in,xmpi_in,xmd_in,xmn_in,xA_in, &
     &  np_in,ne_in,nk_fname_in)
  use mathtool
   implicit none
   integer*8 :: irn_in(nwlk_in)
   integer*4 :: nev_in,nwlk_in,xA_in,i_fg_in,np_in,i,j,ne_in,np0,ne0,ien
   integer*4 :: ipot,i_fsi_in
   real*8 :: xpf_in,xmpi_in,xmd_in,xmn_in,mlept_in,hp,he,thetalept_in
   real*8, allocatable :: pv0(:),dp0(:,:),ep0(:)
   character*40 :: nk_fname_in
   logical :: pwia
   
   nev=nev_in
   nwlk=nwlk_in
   xpf=xpf_in
   xmpi=xmpi_in
   xmd=xmd_in
   xmn=xmn_in
   thetalept=thetalept_in
   xA=xA_in
   i_fg=i_fg_in
   np0=np_in
   ne0=ne_in
   i_fsi=i_fsi_in

   allocate(irn(nwlk))
   irn(:)=irn_in(:)

   if(i_fg.ne.1) then
      open(unit=8,file=nk_fname_in,status='unknown',form='formatted')
      allocate(pv0(np0),dp0(np0,ne0),ep0(ne0))
      do j=1,np0
         read(8,*) pv0(j)!,dp(1,j)         
         read(8,'(4(f6.1,2x,e10.3))')(ep0(i),dp0(j,i),i=1,ne0)
      enddo
      close(8)
   endif
   
   if(i_fg.ne.1) then
      np=2*np0
      ne=ne0
      allocate(pv(np),ep(ne),Pke(np,ne),dp(np))
      hp=pv0(np0)/dble(np)
      ep(:)=ep0(:)
      he=(ep(2)-ep(1))
      do i=1,np
         pv(i)=dble(i-0.5d0)*hp
         do ien=1,ne
            call interpolint(pv0,dp0(:,ien),np0,pv(i),PkE(i,ien),3)
         enddo
      enddo
   else
      np=2*np0
      ne=1
      allocate(pv(np),ep(ne),Pke(np,ne),dp(np))
      hp=xpf/dble(np)
      he=1.0d0
      do i=1,np
         pv(i)=dble(i-0.5d0)*hp
         Pke(i,1)=1.0d0
      enddo
   endif
  
      norm=0.0d0
      do i=1,np
         dp(i)=sum(Pke(i,:))*he
         norm=norm+sum(Pke(i,:))*pv(i)**2*4.0d0*pi*(pv(2)-pv(1))*he
      enddo
      Pke=Pke/norm*4.0d0*pi*xpf**3/3.0d0
      dp=dp/norm*4.0d0*pi*xpf**3/3.0d0
      norm=0.0d0
      do i=1,np
         norm=norm+sum(Pke(i,:))*pv(i)**2*4.0d0*pi*(pv(2)-pv(1))*he
      enddo
      write(6,*) 'norm',norm
      if(pwia)then
         deallocate(Pke,ep)
         allocate(Pke(np,1),ep(1))
          he=1.0d0
          ne=1
          Pke(:,1)=dp(:)
          i_fg=1
      endif
      ! this needs to be updated

      if(i_fsi.eq.1) then
         open(8, file='../pke/realOP_12C_EDAI.dat')
         read(8,*) npot
         allocate(kin(npot),pot(npot))
         do ipot=1,npot
            read(8,*)kin(ipot),pot(ipot)
            pot(ipot)=pot(ipot)
            !....kin and pot are in MeV
         enddo
      endif

      open(10, file='rho_1.dat')
      read(10,*) np_del
      allocate(pdel(np_del),pot_del(np_del))
      do i=1,np_del
         read(10,*) pdel(i),pot_del(i)
      enddo
      
   return
end subroutine   


subroutine mc_eval(ee,w,sig_avg_tot,sig_err_tot)
  use mathtool
  use mympi
   implicit none
   integer*4 :: i,ie1,ie2,ne1,ne2,j
   real*8 :: w,emax,ee
   real*8 :: sig_o(nwlk)
   real*8 :: sig_avg,sig_err
   real*8 :: sig_avg_tot,sig_err_tot

   real*8 :: xpmax,qval_in,sig
   real*8 :: wmax,q2,q2max
   integer*4 :: i_acc,i_avg,i_acc_tot,i_avg_tot
   integer*4 :: ip1_o(nwlk),ie1_o(nwlk),ip2_o(nwlk),ie2_o(nwlk)
   integer*4 :: ip1_n(nwlk),ip2_n(nwlk),ie1_n(nwlk),ie2_n(nwlk)
   real*8 ::  g_o(nwlk),g_n(nwlk),f_o(nwlk),f_n(nwlk) 

   sig_avg=0.0d0
   sig_err=0.0d0
   i_acc=0
   i_avg=0
   g_o=0.0d0
   if (i_fg.eq.1) then
      xpmax=pv(np)
      emax=1.0d0
   else
      xpmax=pv(np)!-pv(1)
      emax=ep(ne)-ep(1)
   endif
   
   !
   do i=1,nwlk
      call setrn(irn(i))
      do while(g_o(i).le.0.0d0)
         ip1_o(i)=1+int(np*ran())
         ie1_o(i)=1+int(ne*ran())
         ip2_o(i)=1+int(np*ran())
         ie2_o(i)=1+int(ne*ran())
         call g_eval(pv(ip1_o(i)),pv(ip2_o(i)),PkE(ip1_o(i),ie1_o(i)), &
              & PkE(ip2_o(i),ie2_o(i)),xpmax,emax,norm,g_o(i))
      enddo
      call getrn(irn(i))
   enddo

   
      
   do i=1,nev
      do j=1,nwlk
         call setrn(irn(j))
         ip1_n(j)=nint(ip1_o(j)+0.05d0*np*(-1.0d0+2.0d0*ran()))
         ie1_n(j)=nint(ie1_o(j)+0.05d0*ne*(-1.0d0+2.0d0*ran()))
         ip2_n(j)=nint(ip2_o(j)+0.05d0*np*(-1.0d0+2.0d0*ran()))
         ie2_n(j)=nint(ie2_o(j)+0.05d0*ne*(-1.0d0+2.0d0*ran()))
         if (ip1_n(j).le.np.and.ip1_n(j).ge.1.and.ie1_n(j).le.nE.and.ie1_n(j).ge.1 &
              &     .and.ip2_n(j).le.np.and.ip2_n(j).ge.1.and.ie2_n(j).le.nE.and.ie2_n(j).ge.1) then
            call g_eval(pv(ip1_n(j)),pv(ip2_n(j)),PkE(ip1_n(j),ie1_n(j)),PkE(ip2_n(j),ie2_n(j)), &
                 &   xpmax,emax,norm,g_n(j))
         else
            g_n(j)=0.0d0
         endif

         if (g_n(j)/g_o(j).ge.ran()) then
            ip1_o(j)=ip1_n(j)
            ip2_o(j)=ip2_n(j)
            ie1_o(j)=ie1_n(j)
            ie2_o(j)=ie2_n(j)
            g_o(j)=g_n(j)
            i_acc=i_acc+1
         endif
         if (i.ge.neq.and.mod(i,nvoid).eq.0) then
            if(ip1_o(j).gt.np.or.ip2_o(j).gt.np) cycle
            call f_eval(ee,pv(ip1_o(j)),pv(ip2_o(j)),ip1_o(j),ip2_o(j),ie1_o(j),ie2_o(j), &
                 & w,sig_o(j))    
            sig_o(j)=sig_o(j)/g_o(j)
            sig_avg=sig_avg+sig_o(j)
            sig_err=sig_err+sig_o(j)**2
            i_avg=i_avg+1
         endif
         call getrn(irn(j))
      enddo
   enddo
   call addall(sig_avg,sig_avg_tot) 
   call addall(sig_err,sig_err_tot)
   call addall(i_avg,i_avg_tot) 
   call addall(i_acc,i_acc_tot) 
   if (myrank().eq.0) then
      sig_avg_tot=sig_avg_tot/dble(i_avg_tot)
      sig_err_tot=sig_err_tot/dble(i_avg_tot)
      sig_err_tot=sqrt((sig_err_tot-sig_avg_tot**2)/dble(i_avg_tot-1))
      write(6,*)'acceptance=',dble(i_acc_tot)/dble(nev*nwlk*nproc())
   endif
 
   return
end subroutine   

subroutine f_eval(ee,p1,p2,ip1,ip2,ie1,ie2,w,sig)
  use mathtool
  implicit none
  integer*4 :: ip1,ip2,ie1,ie2
  real*8 :: ctpp1,p2,ctp2,phip2,p1,ctp1,phip1,ee,eef,w,q2
  real*8 :: qval,cos_theta,jac_c,tan2
  real*8 :: v_ll,v_t
  real*8 :: tnl2,sig0,delta,rho,rhop
  real*8 :: r_now(5),sig,qp(4)

  ctpp1=-1.0d0+2.0d0*ran()
  ctp2=-1.0d0+2.0d0*ran()
  phip2=2.0d0*pi*ran()
  ctp1=-1.0d0+2.0d0*ran()
  phip1=2.0d0*pi*ran()
  sig=0.0d0
  eef=ee/hbarc
  cos_theta=cos(thetalept)
  tan2=(1.0d0-cos_theta)/(1.0d0+cos_theta)
  !.....compute sigma_mott [ fm^2 --> mb ]
  sig0=alpha**2/2.0d0/(1.0d0-cos_theta)/eef**2/tan2
  sig0=sig0*10.0d0
  q2=2.0d0*ee*(ee-w)*(1.0d0-cos_theta)
  if(q2.lt.0.0d0) stop
  qval=sqrt(q2+w**2)
  v_ll=(-q2/qval**2)**2
  v_t=(q2/qval**2/2.0d0+tan2)
  call int_eval(ctpp1,p2,ctp2,phip2,p1,ctp1,phip1,ip1,ip2,ie1,ie2,w,qval,r_now)
  r_now(:)=r_now(:)*2.0d0**3*(2.0d0*pi)**2*ppmax
  sig=sig0*(v_ll*r_now(1)+v_t*r_now(4))*1.e9

  return

end subroutine f_eval

module nuclear_parameters
    ! Initialize all parameters
    ! Dirac matrices
    ! etc.
end module

subroutine one_body(p1, pp1, q, r_now)
    real*8, INTENT IN :: p1(4), pp1(4), q(4)
    real*8, INTENT OUT :: r_now(5)

end subroutine

subroutine two_body(p1, p2, pp1, pp2, q, r_now)
    real*8, INTENT IN :: p1(4), p2(4), pp1(4), pp2(4), q(4)
    real*8, INTENT OUT :: r_now(5)

end subroutine


subroutine int_eval(ctpp1,p2,ctp2,phip2,p1,ctp1,phip1,ip1,ip2,ie1,ie2,w,qval,r_now)
   ! ctpp1: Cos theta_p1'

   ! p2: |p2| momentum
   ! ctp2: Cos theta_p2
   ! phip: Phi_p2

   ! p1: |p1| momentum
   ! ctp1: Cos theta_p1
   ! phip1: Phi_p1

   ! ip1: momentum index of spectral function 
   ! ip2: momentum index of spectral function
   ! ie1: energy index of spectral function 
   ! ie2: energy index of spectral function

   ! w: energy transfer
   ! qval: 3 momentum transfer
   use dirac_matrices         
   use mathtool
   implicit none
   real*8, parameter :: lsq=0.71*1.e6,l3=3.5d0*1.e6,xma2=1.1025d0*1.e6
   real*8, parameter :: fstar=2.13d0,eps=20.0d0
   integer*4 :: ie1,ie2,ip1,ip2
   real*8 :: w,ctpp1,p2,ctp2,phip2,p1,ctp1,phip1,stpp1,stp1,stp2
   real*8 :: at,bt,vt,par1,par2,pp1,den,jac,arg,qval
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
end subroutine   




subroutine g_eval(pj1,pj2,pke1,pke2,pmax,emax,norm,g)
    implicit none
    real*8, parameter :: pi=acos(-1.0d0)
    real*8 :: pj1,pj2,pke1,pke2,g,pmax,emax,norm
    g=(4.0d0*pi)*pj1**2*pke1
    g=g*(4.0d0*pi)*pj2**2*pke2
    g=g/norm**2
!    g=g*(1.0d0+cost/2.0d0)/2.0d0
!    g=g*(1.0d0+tmu/tmax)/(3.0d0/2.0d0*tmax)
    
    return
    end subroutine


end module        
