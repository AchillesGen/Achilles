    module xsec_fact
    implicit none
    real*8, parameter :: xmp=938.272046d0,xmn=939.56563d0,xm=0.5d0*(xmn+xmp)
    real*8, parameter :: hbarc=197.327053d0,G_F = 1.14e-5*0.1973d0**2
    real*8, parameter :: pi= 4.0d0*atan(1.0d0),c0=1.0d0/16.0d0/pi/pi
    real*8, save :: e,ef,cost,q2,p2,pf2,sint,sina2,wf,wtf
  end module xsec_fact



  subroutine cc1(xq,w,wt,xk,xp,phi,ee0,theta,ig,sig)
    use xsec_fact
    use dirac_matrices
    implicit none
    integer*4 :: ig,i
    real*8 :: xq, w,wt, xk, xp, ee0, theta,sig
    real*8 :: ek,epf,xmf
    real*8 :: phi,xk_x,xk_y,xk_z
    real*8 :: p_4(4),pp_4(4),q_4(4), cosa
!     ee0      incident electron energy in MeV
!     w        electron energy loss in MeV
!     xq       three-momentum transfer in inverse fm
!     thetae   electron scattering angle
!     xk, xp   initital and final nucleon three-momenta in inverse fm
!     sina2    sin(alpha)**2, alpha being the angle between the 
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
      use xsec_fact
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






  

