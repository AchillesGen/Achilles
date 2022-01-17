    module xsec_fact
    implicit none
    real*8, parameter :: xmp=938.272046d0,xmn=939.56563d0,xm=0.5d0*(xmn+xmp)
    real*8, parameter :: hbarc=197.327053d0,G_F = 1.14e-5*0.1973d0**2
    real*8, parameter :: pi= 4.0d0*atan(1.0d0),c0=1.0d0/16.0d0/pi/pi
    real*8, save :: e,ef,cost,q2,p2,pf2,sint,sina2,wf,wtf
  end module xsec_fact



  subroutine cc1(xq,w,wt,xk,xp,p_4,pp_4,k_4,kp_4,theta,ig,sig_p, sig_n)
    use xsec_fact
    use dirac_matrices
    implicit none
    integer*4 :: ig,i
    real*8 :: xq, w,wt, xk, xp, ee0, theta
    real*8, intent(out) :: sig_p, sig_n
    real*8 :: ek,epf,xmf
    real*8 :: p_4(4),pp_4(4),q_4(4),cosa
    real*8 :: k_4(4), kp_4(4),kdotq
    
!     ee0      incident electron energy in MeV
!     w        electron energy loss in MeV
!     xq       three-momentum transfer in inverse fm
!     thetae   electron scattering angle
!     xk, xp   initital and final nucleon three-momenta in inverse fm
!     sina2    sin(alpha)**2, alpha being the angle between the 
      wf = w/hbarc
      wtf=wt/hbarc
      ef = kp_4(1)           
      q2 = xq**2
      q_4(1)=wtf
      q_4(2:3)=0.0d0
      q_4(4)=xq
      !
      !write(6,*) 'initial lepton',k_4(1:4)
      !write(6,*) 'final lepton',kp_4(1:4)
      
      call current_init(p_4,pp_4,q_4,k_4,kp_4,wf)
      
      call define_spinors()
      call define_lept_spinors()
      
      call sigccc(1,sig_p,ig)
      call sigccc(2,sig_n,ig)
      sig_p=sig_p*ef**2
      sig_n=sig_n*ef**2
      
      return
    end subroutine cc1
    
    subroutine sigccc(id,sig,ig)
      use xsec_fact
      use dirac_matrices
      implicit none
      integer*4 :: ig,id
      real*8, parameter :: alpha=1.0d0/137.0d0
      real*8 :: tan2,sig_mott,qm2
      real*8 :: ff1s,ff2s,ff1v,ff2v,ffa,ffp,ges,gms,gev,gmv
      real*8 :: sig,sig_p,sig_n,ff1p,ff1n,ff2p,ff2n


      qm2 = q2-wf**2           
      sig_mott = alpha**2/qm2**2*10.0d0*4.0d0/2.0d0

      call nform(ig,qm2,ff1s,ff2s,ff1v,ff2v,ffa,ffp,ges,gms,gev,gmv)
      ff1p=0.5d0*(ff1v+ff1s)
      ff2p=0.5d0*(ff2v+ff2s)
      ff1n=0.5d0*(-ff1v+ff1s)
      ff2n=0.5d0*(-ff2v+ff2s)

      if(id.eq.1) then
         call det_Ja(ff1p,ff2p)
         call contract(sig_p)
         sig=sig_mott*0.5d0*(sig_p) 
      elseif(id.eq.2) then
         call det_Ja(ff1n,ff2n)
         call contract(sig_n)
         sig=sig_mott*0.5d0*(sig_n) 
      endif
      return
    end subroutine sigccc






  


