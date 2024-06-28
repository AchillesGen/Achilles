! Inputs
!
! four-momenta in MeV in the original frame (Lab, CM, etc.)
!       qvec_in                   ! virtual photon four-momentum (thie needs to be along z-axis)
!       xpnuc_in                  ! final nucleon four-momentum
!       xk_in                     ! final meson four-momentum

!       mode_in =1  neutrino C! scattering on nucleon
!               =3  anti-neutrino C! scattering on nucleon
!               =-1  neutrino N! scattering on nucleon
!               =-3  anti-neutrino N! scattering on nucleon
!               =10 electron scattering on nucleon
!       isl                     ! =0: DCC; =1:SL
!       itiz                    ! incoming nucleon isospin z-component * 2; =1:p ; =-1:n
!       itpiz                   ! isospin z-component for final meson
!       ipar                    ! =1 :pion ; =2: eta
!
!
! Output
!       zj_mu(isf,isi,ig)         ! current matrix elements in the original frame
!             isf                 ! 2* (z-component of final nucleon spin)
!             isi                 ! 2* (z-component of incoming nucleon spin)
!             ig                  ! 0,1,2,3 : photon polarization
! Relation between zj_mu and 
! gamma N -> piN differentianl cross section (d\sigma / d\Omega_\pi [\mu b/sr]) in CM
! when qvec_in, xpnuc_in, and xk_in are given in pi-nucleon CM frame
!       dsigma : d\sigma / d\Omega_\pi [\mu b/sr] 
!       wcm    : piN invariant mass
!       fnu!   : nucleon mass
!       xkabs  : magnitude of pion momentum (spatial component)
!       alfa   : the fine structure constant  1/137
!       pi     : 3.1415...
!       E_gamma=(wcm**2-fnuc**2)/(2*fnuc)
!       dsigma=0
!       do isi=-1,1,2
!       do isf=-1,1,2
!       do ig=1,2
!       dsigma=dsigma
!      &  + 4*pi**2*alfa/E_gamma
!      &    *abs(zj_mu(isf,isi,ig))**2
!      &    *fnuc*xkabs/(16*pi**3*wcm)
!      &    /4                    ! initial spins averaged
!      &    *1d4                  ! \mu b/sr
!       end do
!       end do
!       end do




      subroutine set_param
      implicit real*8 (a-h,o-y)
      implicit complex*16 (z) 
      parameter (maxmb=1)
      common / cscal / pi,fm
      common/masses/am1dt(maxmb),am2dt(maxmb)
      common /icgmb / icgmb(5)
      !fm      = 197.3289d0
      !pi = atan(1d0)*4d0
      !alfa = 1/137.0359895d0

      imwpon=0
      isw0=0
      isigma=1

      icgmb    = 0
      icgmb(1) = 1 ! igmb=1 pn
      icgmb(2) = 2 ! igmb=2 en
      icgmb(3) = 6 ! igmb=3 on
      icgmb(4) = 7 ! igmb=4 kl
      icgmb(5) = 8 ! igmb=5 ks
!=============================================
!    DEFINE MASSES FOR CHANNELS
!=============================================

!     1 :pi-N      2 :eta-N      3 :pi-Delta
!     4 :sigma-N   5 :rho-N      6 :omega-N
!     7 :K-Lambda  8 :K-Sigma
!
      amn   = 938.5d0
      api   = 138.5d0
!      aeta  = 547.45d0
      aeta  = 548.d0
      amdel = 1299.d0
!      asigma= 896.8d0
!      arho  = 811.7d0
      asigma= 897.d0
      arho  = 812.d0

      aomega= 782.d0  
      akaon = 495.d0  
      amlam = 1115.7d0 
      amsig = 1193.d0 

      if(isw0.eq.1) then   !JLMS
      amn   = 940.d0
      api   = 140.d0
      aeta  = 550.d0
      amdel = 1299.d0
      asigma= 897.d0
      arho  = 812.d0
      end if
      if(isigma.eq.1) then !new D->piN,sig->pipi f.f.
      amdel = 1280.d0
      asigma= 700.d0
      end if

!mwp's values
      if(imwpon.eq.1) then
      amn   = 938.5d0
      api   = 138.5d0
      aeta  = 547.5d0
      amdel = 1300.d0
      asigma= 898.6d0
      arho  = 811.7d0
      aomega= 782.6d0
      end if

      am1dt(1) = api
      am2dt(1) = amn
      do i=1,1
      am1dt(i) = am1dt(i)/fm
      am2dt(i) = am2dt(i)/fm
      end do

      return
      end




      subroutine pion_init(fneu,fpro,feta_in,fpio0,fpi1,pi_in,fm_in)
       implicit real*8(a-h,o-y)
       common / ccoup / fpio,fnuc,flep,fnuci,fdeu,feta,fmfin
       common / cscal / pi,fm


       !fneu = 939.56563d0/fm
       !fpro = 938.27231d0/fm
       !feta =  548.d0/fm
       !fpio1 = 139.57018d0/fm
       !fpio0 = 134.9764d0/fm
       fnuc=(fpro+fneu)/2
       !fpio=(2*fpio1+fpio0)/3
       fpio = fpio0
       fnuc2=fnuc**2
       fpio2=fpio**2
       fnuci = fnuc
       fnuci2=fnuci**2
       feta=feta_in
!...
       fm = fm_in!197.3289d0
       pi = pi_in!atan(1d0)*4d0


       write(6,*)'fneu= ', fneu*fm 
       write(6,*)'fpro = ', fpro*fm 
       write(6,*)'fnuc = ', fnuc*fm 
       write(*,*)'fpio = ', fpio*fm

      end subroutine
      
       




      subroutine amplitude(qvec_in,xpnuc_in,xk_in,mode_in,isl,itiz,itpiz &
     &  ,ipar,zj_mu)
      implicit real*8(a-h,o-y)
      implicit complex*16(z)
      common / cscal / pi,fm
      common / ccoup / fpio,fnuc,flep,fnuci,fdeu,feta,fmfin
      parameter(ndimd=300)
      parameter (maxwcm=70,maxq2=30,npar=2)
      common/ rfac / rfac 
      common/ mode / vfac,vvfac(-1:1),mode
      common / igmbs / tm_f,igmbs(npar),igmb_final_meson
      parameter (maxmb=1 ,maxlsj=20,maxmom=1)
      common/masses/am1dt(maxmb),am2dt(maxmb)
      common /icgmb / icgmb(5)
      dimension isbi(8),ismi(8),ismix(8)
      data isbi/1,-1,1,-1,1,-1,1,-1/
      data ismi/1,1,0,0,-1,-1,2,2/
      data ismix/1,1,0,0,-1,-1,0,0/
      parameter (njmx=5,lpimax=5)
      dimension dfun(2*njmx-1,-2*njmx+1:2*njmx-1,-2*njmx+1:2*njmx-1)
      dimension zphi_pins(-lpimax:lpimax),zphi_qs(-2*njmx-2:2*njmx+2)
     &  ,fn_sign(-2:2)
      common / zmtx / zmtx(8,maxlsj,5,5)
      common/chdat1/njLs,jpind(maxlsj) ,Lpind(maxlsj)
     &                  ,ispind(maxlsj),itpind(maxlsj)
      dimension zjx_mu(-1:1,-1:1,0:3),zj_mu(-1:1,-1:1,0:3)
      common / nLsdt / nLsdt(maxmb,maxlsj),Ldt0(5,maxmb,maxlsj)
      parameter(nchn=5,npmx=30)
      dimension icnv(maxmb),xk2cm(0:3),xlrs(0:3,0:3),pnuc2cm(0:3)
      data icnv/1/
      dimension pcm(0:3),xlr(0:3,0:3)
     &  ,xpnuc_in(0:3),xk_in(0:3),qvec_in(0:3)
     &  ,xpnuc(0:3),xk(0:3),qvec(0:3)
     &  ,bleg_pin(0:8,-8:8),itizs(2),sym(2)
     &  ,pnucx(0:3)
      dimension  px(0:3),ixi_cnv(8),px2cm(0:3)
     &  ,zmxpi(2,2),zmxpf(2,2),qx2cm(0:3)
      dimension isp(-1:1),zcrnt(2,2,-1:2)
     &  ,zcrntx(2,2,-1:2),zsumx(2,2)
      common / eps_tpin / eps_tpin
      data init /1/
      data ixi_cnv /1,2,5,6,3,4,7,8/
      data isp /2,0,1/
      save init,zei,sqhf,lmax,tpiz,fnuc2,fnuci2
     & ,facp,tcrz,fac
      save tmax,imxi,igm1_max,igm1_step,fn_sign
      save irot_spin,irot_q,irotq_sign,icomp,num_pol1,num_pol2


c
      if(init==1)then
      init=0
      zei=(0,1d0)
      sqhf=1/sqrt(2d0)
      lmax=5
      imxi=8
      igm1_max=2
      igm1_step=1
      icomp=0
c
      mode=mode_in
      imxi=8
      igm1_max=2
      igm1_step=1
      num_pol1=0
      num_pol2=3
!
      irot_spin=1               ! =1 : spin rotation ; =0: no spin rotation
      irot_q=1     ! =0 photon momentum in 2CM is along z-axis (approx); =1 not along z-axis
      
      irotq_sign=(-1)**(irot_q+1)
!
      do i1=-2,2,2
      fn_sign(i1)=(-1)**(i1/2)
      end do
!
      call set_param
      call bifc
!      
      rfac=1
c
      if(ipar==1)then
        fmfin=fpio
        igmb_final_meson=1
        tm_f=1d0
      else if(ipar==2)then
        fmfin=feta
        igmb_final_meson=2
        tm_f=0d0
        itpiz=0
      else
        stop 'ipar error'
      endif  
!
      call read_amp(isl,icomp,ipar,0,0  ! called after fnuci is set
     & )
!        
      tmax=tm_f+.5d0+eps_tpin
!
      if(isl==1)then
      facp=1d0/fnuc2
      else
      facp=fm**2
      endif 
      fac=facp
     &  /(4*pi)
      fac=sqrt(fac)
      if(mode.ge.1.and.mode.le.4)fac=fac*sqrt(2d0) ! isospin factor for CC
      endif                     ! init=1
      tpiz=itpiz
      if(abs(itpiz).gt.1)stop 'itpiz error'
!
      if(mode==1.or.mode==2)then
        tcrz=1
      elseif (mode==3.or.mode==4)then
        tcrz=-1
      elseif (mode.ge.-4.and.mode.le.-1.or.mode==10.or.mode==11)then
        tcrz=0
      endif
!
      if(mode.ge.1)then         ! C! and EM
        vfac = 1
      else if(mode.le.-1)then   ! NC
        sw2 = 0.2312d0          ! Weinberg angle:  sin^2 \theta_W
        vfac = 1 - 2*sw2
        vvfac(1) =  - 2*sw2
        vvfac(-1) =    2*sw2
      endif 
!
      xpnuc=xpnuc_in/fm
      xk=xk_in/fm
      qvec=qvec_in/fm
!
      wcm2=(xk(0)+xpnuc(0))**2-(xk(1)+xpnuc(1))**2
     &  -(xk(2)+xpnuc(2))**2-(xk(3)+xpnuc(3))**2
      Q2=qvec(1)**2+qvec(2)**2+qvec(3)**2-qvec(0)**2
!
! impulse
!
      zj_mu=0
      tiz=dble(itiz)/2
      tpinz=tcrz+tiz
      if(abs(tpinz-tpiz).gt.0.50001)
     &  stop 'isospin (input parameter) mismatch'
!
      do ix=0,3
      pcm(ix)=xpnuc(ix)+xk(ix)
      pnucx(ix)=xpnuc(ix)
      px(ix)=xpnuc(ix)+xk(ix)-qvec(ix) ! incoming nucleon momentum 
      end do
      px(0)=sqrt(fnuc2+px(1)**2+px(2)**2+px(3)**2) ! on mass shell is assumed
      call lorentz_trans(1,1,pcm,xlrs) ! Lab -> 2CM
!
      do ix1=0,3
      xk2cm(ix1)=0
      px2cm(ix1)=0
      qx2cm(ix1)=0
      do ix2=0,3
      xk2cm(ix1)=xk2cm(ix1)+xlrs(ix1,ix2)*xk(ix2)
      px2cm(ix1)=px2cm(ix1)+xlrs(ix1,ix2)*px(ix2)
      qx2cm(ix1)=qx2cm(ix1)+xlrs(ix1,ix2)*qvec(ix2)
      end do
      pnuc2cm(ix1)=-xk2cm(ix1)
      end do
      pnuc2cm(0)=sqrt(fnuc2+pnuc2cm(1)**2+pnuc2cm(2)**2+pnuc2cm(3)**2)
      xk2cmabs=sqrt(xk2cm(1)**2+xk2cm(2)**2+xk2cm(3)**2)
      qx2cmabs=sqrt(qx2cm(1)**2+qx2cm(2)**2+qx2cm(3)**2)
!
      xz_pin=xk2cm(3)/xk2cmabs
      sxz_pin=sqrt(1-xz_pin**2)
      cphi_pin=xk2cm(1)/(sxz_pin*xk2cmabs)
      sphi_pin=xk2cm(2)/(sxz_pin*xk2cmabs)
      zphi_pin=cphi_pin+zei*sphi_pin
      zphi_pins(0)=1
      zphi_qs(0)=1
      if(irot_q==0)then
      do llz=1,2
      zphi_pins(llz)=zphi_pin**llz
      zphi_pins(-llz)=1d0/zphi_pins(llz)
      end do
      else                      ! irot_q=1
      do llz=1,lpimax
      zphi_pins(llz)=zphi_pin**llz
      zphi_pins(-llz)=1d0/zphi_pins(llz)
      end do
      xz_q=qx2cm(3)/qx2cmabs
      sxz_q2=1-xz_q**2
      if(sxz_q2.gt.0)then
      sxz_q=sqrt(sxz_q2)
      cphi_q=qx2cm(1)/(sxz_q*qx2cmabs)
      sphi_q=qx2cm(2)/(sxz_q*qx2cmabs)
      zphi_q=cphi_q+zei*sphi_q
      else
      xz_q=1d0
      sxz_q=0
      zphi_q=(1d0,0d0)
      endif 
      do mj=1,njmx+1
      zphi_qs(2*mj)=zphi_q**mj
      zphi_qs(-2*mj)=1d0/zphi_qs(2*mj)
      end do
      call setdfun(xz_q,dfun)
      endif                     ! irot_q
!
      call ylmsub(lmax,xz_pin,bleg_pin)
      call lorentz_trans(2,1,pcm,xlr) ! 2CM -> Lab
!
      wcm=sqrt(wcm2)
      call interpolate_amp(wcm,Q2,icomp,itiz,tpinz)
!
      do 131 igmb = igmb_final_meson,igmb_final_meson
      if(icgmb(igmb).eq.0) go to 131
      ic        = icgmb(igmb)
      am1=am1dt(ic)
      am2=am2dt(ic)
      if(wcm.lt.am1+am2)goto 131

      zcrnt=0
      do ixi1p = 1,imxi
      ixi1 = ixi_cnv(ixi1p)
      igm1 = ismi(ixi1)
      igm1x = ismix(ixi1)
      lambda_N=-isbi(ixi1)
      lamx=isp(lambda_N*irotq_sign)
      Lambda_i=2*igm1x-lambda_N
      do i=1,njLs
      jpin=jpind(i)
      Lpin=Lpind(i)
      itpin=itpind(i)
      tpin=dble(itpin)/2
      xlpin=dble(Lpin)/2
      xjpin=dble(jpin)/2
      llpin=Lpin/2
      if(tpin+eps_tpin.gt.abs(tpinz).and.tmax.gt.tpin
     &    .and.jpin.ge.abs(Lambda_i))then
      do ils=1,nLsdt(ic,i)
      zfac_a=sqrt(dble(jpin+1))
     &    *cbg(1d0,tcrz,0.5d0,tiz,tpin,tpinz)
     &    *cbg(tm_f,tpiz,0.5d0,tpinz-tpiz,tpin,tpinz)
     &    *zmtx(ixi1,i,igmb,ils)
      do isf=-1,1,2
      isfx=isp(isf)
      xs=dble(isf)/2
      if(irot_q==0)then
      mj=Lambda_i
      xmj=dble(mj)/2
      llz=(mj-isf)/2
      zzz=cbg(xlpin,xmj-xs,0.5d0,xs,xjpin,xmj)
     &    *bleg_pin(llpin,llz)
     &    *zphi_pins(llz)
      else 
      zzz=0
      do mj=max(-Lpin+isf,-jpin),min(Lpin+isf,jpin),2
      xmj=dble(mj)/2
      llz=(mj-isf)/2
      zzz=zzz
     &    +cbg(xlpin,xmj-xs,0.5d0,xs,xjpin,xmj)
     &    *bleg_pin(llpin,llz)
     &    *zphi_pins(llz)
     &    *dfun(jpin,mj,Lambda_i)
     &    *zphi_qs(-mj+Lambda_i)
      end do                    ! mj
      endif 
      zcrnt(isfx,lamx,igm1)=zcrnt(isfx,lamx,igm1)
     &    +zfac_a*zzz
      end do                    ! isf
      end do                    ! ils
      endif                     ! tpin.ge.abs(tpinz)
      end do                    ! i
      end do                    ! ixi1p
! Conversion from helicity basis to spin basis
      if(irot_q==1)then
      zcrntx=zcrnt
      do igm1=-1,igm1_max,igm1_step
      do isi=-1,1,2
      isix=isp(isi)
      do isfx=1,2
      zcrnt(isfx,isix,igm1)=0
      do lambda_N=-1,1,2
      lamx=isp(lambda_N)
      zcrnt(isfx,isix,igm1)=zcrnt(isfx,isix,igm1)
     &  +dfun(1,-lambda_N,isi)
     &  *fn_sign(lambda_N+isi)*zphi_qs(isi+lambda_N)
cc!     &    *(-1)**((1-isi)/2)*zphi_qs(isi-lambda_N) ! modified as above to account for the used base
     &    * zcrntx(isfx,lamx,igm1)
      end do
      end do
      end do
      end do
      endif                     ! irot_q=1
! spin rotation
      if(irot_spin==1)then
      zcrntx=zcrnt
!     for the D-function (zmxpi) of the incoming gamma-N
      call rspin(pcm,px,px2cm,zmxpi)
!     for the D-function(zmxpf) of the final pi-N
      call rspin(pcm,pnucx,pnuc2cm,zmxpf)
      zcrnt=0
      do igm1=-1,igm1_max,igm1_step
      zsumx=0
      do i4=1,2
      do i2=1,2
      do i3=1,2
      zsumx(i3,i2)=zsumx(i3,i2)+zcrntx(i3,i4,igm1)*zmxpi(i2,i4)
      end do
      end do
      end do
      do i3=1,2
      do i2=1,2
      do i1=1,2
      zcrnt(i1,i2,igm1)=zcrnt(i1,i2,igm1)
     &    +conjg(zmxpf(i1,i3))*zsumx(i3,i2)
      end do
      end do
      end do
      end do                    ! igm1
      endif                     ! irot_spin=1
  131 continue                  ! igmb

      do is1=-1,1,2
      do is2=-1,1,2
      is1x=isp(is1)
      is2x=isp(is2)
      zjx_mu(is1,is2,0)=zcrnt(is1x,is2x,0)
      zjx_mu(is1,is2,3)=zcrnt(is1x,is2x,2)
      zjx_mu(is1,is2,1)
     &  =(zcrnt(is1x,is2x,-1)-zcrnt(is1x,is2x,1))*sqhf
      zjx_mu(is1,is2,2)
     &  =(zcrnt(is1x,is2x,-1)+zcrnt(is1x,is2x,1))*sqhf*zei
      do ix1=0,3
      do ix2=num_pol1,num_pol2
      zj_mu(is1,is2,ix1)=zj_mu(is1,is2,ix1)
     &    +xlr(ix1,ix2)*zjx_mu(is1,is2,ix2)
      end do
      end do
      end do                    ! is2
      end do                    ! is1
      zj_mu=zj_mu*fac
      return
      end

      

      subroutine read_amp(isl,icomp,ipar,nnresc,npiresc
     &  )
      implicit real*8 (a-h,o-y)
      implicit complex*16 (z)
      parameter (maxmb=1 ,maxlsj=20,ndimd=300,maxmom=1)
      common/masses/am1dt(maxmb),am2dt(maxmb)
      common / cscal / pi,fm
      common /icgmb / icgmb(5)
      common / zmtx / zmtx(8,maxlsj,5,5)
      common/chdat1/njLs,jpind(maxlsj) ,Lpind(maxlsj)
     &                  ,ispind(maxlsj),itpind(maxlsj)
      parameter (maxwcm=70,maxq2=30,ilmax=1,mbs=2)
      common / w_q2 /wcms(maxwcm),q2s(maxq2),maxw,mxq2
      common / amp_offshell /
     &   zampv(maxq2,maxwcm,maxmom,8,maxlsj,mbs,ilmax)
     &  ,zampv_is(maxq2,maxwcm,maxmom,8,maxlsj,mbs,ilmax)
     &  ,zampa(maxq2,maxwcm,maxmom,8,maxlsj,mbs,ilmax)
      parameter(ndim=60)
      dimension xin(ndim),yin(ndim),xout(ndim),yout(ndim)
      dimension xin2(ndim),yin2(ndim),yout2(ndim),xout2(ndim)
      dimension yin3(ndim),yin4(ndim)
      common / nLsdt / nLsdt(maxmb,maxlsj),Ldt0(5,maxmb,maxlsj)
      dimension is1dt(maxmb),is2dt(maxmb)
      data is1dt/0/,is2dt/1/
      common / ccoup / fpio,fnuc,flep,fnuci,fdeu,feta,fmfin
      parameter(nchn=5,npmx=30)
      common /para_dfn / jmax
      parameter (njmx=5,lpimax=5)
      common / cdfi  / meshx,mxx,mxj,mxm
      dimension za(4),idata(3)
      common/ i_pipin/ i_pipin,ic_ppn1,ic_ppn2
      common/ rfac / rfac      
      common/ mode / vfac,vvfac(-1:1),mode
      dimension icnv(maxmb)
      data icnv/1/
      parameter (ndeu=500,npar=2)
      dimension xinx(ndeu),yinx(ndeu),yinx2(ndeu),xoutx(ndeu)
     &  ,youtx(ndeu),youtx2(ndeu)
      common/ zamp_pin /zamp_pin(maxwcm,ndeu,npar,maxlsj),pmom(ndeu),mxp
      common / zamp_pin_save / zamp_pin_save(maxmom,npar,maxlsj)
      common / eps_tpin / eps_tpin
      common / igmbs / tm_f,igmbs(npar),igmb_final_meson
      dimension id1(4),id2(4),itizs(2),cgiso(-1:1,3)
     &  ,pmom_mid(maxmom)
      data id1 /1,2,3,7/
      data id2 /6,5,4,8/
      save modex,pmom_mid,i_offshell2,tmax,ndir
      save wcm_old,Q2_old,itiz_old,tpinz_old

      modex=1
      eps_tpin=0.01
      tmax=tm_f+.5d0+eps_tpin

      if(isl==1)then
      open(10,file='sl_EW.dat',form='formatted',status='old')
      else
      open(10,file="data/dcc_EW.dat",form='formatted',status='old')
      endif 
      read(10,*) njLs
      jmax=0
      lmax=0
      do i=1,njLs
      read(10,*) jpind(i) ,Lpind(i),ispind(i),itpind(i)
      jmax=max(jmax,jpind(i))
      lmax=max(lmax,Lpind(i))
      end do

      if(jmax.gt.2*njmx-1)stop 'njmx error'
      if(lmax.gt.2*lpimax)stop 'lpimax error'
!      
! set nLsdt
      do i  =1,njLs
      jpin  =jpind(i)
      Lpin  =Lpind(i)
      ispin =ispind(i)
      ipar00=(-1)**(Lpin/2 +2+1)

      do 50 ic=1,maxmb
      is1=is1dt(ic)
      is2=is2dt(ic)
      ismin=abs(is1-is2)
      ismax=is1+is2
      ni=0
      do 51 is=ismin,ismax,2
      Lmin=abs(jpin-is)
      Lmax=jpin+is
      do 51 L=Lmin,Lmax,2
      ipa=(-1)**(L/2+2+1)
      if(ic.eq.4)ipa=(-1)**(L/2+2)
      if(ipa.eq.ipar00)then
      ni=ni+1
      Ldt0(ni,ic,i) = L
      end if
   51 continue
      nLsdt(ic,i)  = ni
   50 continue       
      end do                    ! end loop i

      read(10,*) maxw,mxq2
      do i=1,maxw
      read(10,*) wcms(i)
      wcms(i)=wcms(i)/fm
      end do
      do i=1,mxq2
      read(10,*) q2s(i)
      q2s(i)=q2s(i)/fm**2        
      end do
      if(maxw.gt.ndeu)stop 'ndeu error 2'
      if(maxw.gt.maxwcm)stop 'maxwcm error'
      if(mxq2.gt.maxq2)stop 'maxq2 error'

      idata(1)=1                ! bare
      idata(2)=1                ! dressed
      idata(3)=1                ! non-res      
! stable final states
      mxp=0
      ip=mxp+1
      read(10,*) namp1,namp2,namp3
      do ix=1,namp1
      read(10,2800,end=100)ie,iq,idx,ipw,igmb,ils
     &    ,za(1),za(2),za(3)
! za(i) i=1:bare; 2:dressed N*; 3: non-res
      if(igmb.le.mbs.and.ipw.le.njLs)then
      zampv(iq,ie,ip,idx,ipw,igmb,ils)=0
      do n=1,3
      zampv(iq,ie,ip,idx,ipw,igmb,ils)
     &    =zampv(iq,ie,ip,idx,ipw,igmb,ils)+za(n)*idata(n)
      end do
      endif 
      end do
  100 continue

      isign=1
      if(mode==10.or.mode==11)isign=-1 ! correct phase for neutron amp; isospin CG multiplied later
      do ix=1,namp2
        read(10,2800,end=106)ie,iq,idx,ipw,igmb,ils
     &    ,za(1),za(2),za(3)
      if(igmb.le.mbs.and.ipw.le.njLs)then
      zampv_is(iq,ie,ip,idx,ipw,igmb,ils)=0
      do n=1,3
      zampv_is(iq,ie,ip,idx,ipw,igmb,ils)
     &    =zampv_is(iq,ie,ip,idx,ipw,igmb,ils)+za(n)*idata(n)*isign
      end do
      endif 
      end do
  106 continue

      if(mode.lt.10)then        ! neutrino case
      do ix=1,namp3
      read(10,2800,end=110)ie,iq,idx,ipw,igmb,ils
     &    ,za(1),za(2),za(3)
      if(igmb.le.mbs.and.ipw.le.njLs)then
      zampa(iq,ie,ip,idx,ipw,igmb,ils)=0
      do n=1,3
      zampa(iq,ie,ip,idx,ipw,igmb,ils)
     &    =zampa(iq,ie,ip,idx,ipw,igmb,ils)+za(n)*idata(n)
      end do
      endif 
      end do
  110 continue 
      endif 
 2800 format (6i4,10g20.10)
      close(10) 

!conversion 1/2p 1/2n -> 1/2v 1/2s basis

      if(mode.lt.10)then        ! neutrino case
      do 510 ipw=1,njLs
      if(itpind(ipw)==3)goto 510
      do 500 ie=1,maxw
      do 500 iq=1,mxq2
      do 500 idx=1,3
      do 500 igmb=1,mbs
      ic        = icgmb(igmb)
      do 500 ils=1,nLsdt(ic,ipw)
      zp=(zampv(iq,ie,ip,idx,ipw,igmb,ils)
     &    +zampv_is(iq,ie,ip,idx,ipw,igmb,ils))*0.5d0
      zm=(zampv(iq,ie,ip,idx,ipw,igmb,ils)
     &  -zampv_is(iq,ie,ip,idx,ipw,igmb,ils))*0.5d0
      zampv(iq,ie,ip,idx,ipw,igmb,ils)=zm    ! isovector
      zampv_is(iq,ie,ip,idx,ipw,igmb,ils)=zp ! isoscalar
  500 continue 
  510 continue
      endif                     ! mode < 10

! z-component from time component via current conservation for vector current
! V_z = V_0 * omega_cm / q_cm
      do 535 ie=1,maxw
        wcmx=wcms(ie)
      do 535 iq=1,mxq2
        Q2x=q2s(iq)
        qc0=(wcmx**2-fnuci**2-Q2x)/(2*wcmx)
        qc2=Q2x+qc0**2
        qc=sqrt(qc2)
        xxx=qc0/qc
      do 530 ipw=1,njLs
      do 530 idx=3,4
        idxx=idx+4
      do 530 igmb=1,mbs
      ic        = icgmb(igmb)
      do 530 ils=1,nLsdt(ic,ipw)
      zampv(iq,ie,ip,idxx,ipw,igmb,ils)
     &    =zampv(iq,ie,ip,idx,ipw,igmb,ils)*xxx
      zampv_is(iq,ie,ip,idxx,ipw,igmb,ils)
     &  =zampv_is(iq,ie,ip,idx,ipw,igmb,ils)*xxx
  530 continue 
  535 continue 

      return

! note : interpolation of structure function to be done

! production amp for IA diagrams
      entry interpolate_amp(wcm,Q2,icomp,itiz,tpinz)

      if(abs(wcm-wcm_old).lt.1d-5.and.abs(Q2-Q2_old).lt.1d-5
     &  .and.itiz==itiz_old.and.abs(tpinz-tpinz_old).lt.1d-5)return
      wcm_old=wcm
      Q2_old=Q2
      itiz_old=itiz
      tpinz_old=tpinz

      iax=1                     ! axial
      ivec=1                    ! vector
      i_pion_pole=1             ! pion-pole term in axial current
      zei=(0,1d0)
      ip=mxp+1
      idxp_mx=4
      idxp_mx_v=3

      qc0=(wcm**2-fnuci**2-Q2)/(2*wcm)
      qc2=Q2+qc0**2
      qc=sqrt(qc2)
      xxx=qc0/qc

      if(Q2==0d0)Q2=1d-6/fm**2
      do i=1,maxw
        if(wcm.le.wcms(i))goto 200
      end do
      write(*,*) wcm*fm,wcms(maxw)*fm
      stop 'wcm interpolation error'
  200 iw=i
      do i=1,mxq2
        if(Q2.le.q2s(i))goto 210
      end do
      write(*,*) Q2*fm**2*1d-6,q2s(mxq2)*fm**2*1d-6
      stop 'Q2 interpolation error'
  210 iq=i
      if(iw==1)then
        write(*,*) wcm*fm
        stop 'wcm interpolation error 2'
      endif 
      if(iq==1)then
        write(*,*) Q2*fm**2*1d-6
        stop 'Q2 interpolation error 2'
      endif

      if(iw.le.2)then
        iw_start=1
        iw_end=4
      else if(iw.ge.maxw-1)then
        iw_start=maxw-3
        iw_end=maxw
      else
        iw_start=iw-2
        iw_end=iw+1
      endif 
      if(iq.le.2)then
        iq_start=1
        iq_end=4
      else if(iq.ge.mxq2-1)then
        iq_start=mxq2-3
        iq_end=mxq2
      else
        iq_start=iq-2
        iq_end=iq+1
      endif 

      iwx=0
      do iw=iw_start,iw_end
      iwx=iwx+1
      xin(iwx)=wcms(iw)
      end do
      iqx=0
      do iq=iq_start,iq_end
      iqx=iqx+1
      xin2(iqx)=q2s(iq)
      end do
      xout(1)=wcm
      xout2(1)=q2

      if(mode.lt.10)then
! axial current interpolation
      do 310 ipw=1,njLs
      tpin=dble(itpind(ipw))/2
      if(tpin+eps_tpin.gt.abs(tpinz).and.tmax.gt.tpin)then
        phv=(-1)**((jpind(ipw)-1)/2 + Lpind(ipw)/2+1)
        pha=-phv
        if(jpind(ipw)==1)then
          idxp_start=2
        else
          idxp_start=1  
        endif 
      do 300 idxp=idxp_start,4
        idx=id1(idxp)
        idxx=id2(idxp)
      do 300 igmb=igmb_final_meson,igmb_final_meson
      ic        = icgmb(igmb)
      do 300 ils=1,nLsdt(ic,ipw)
! stable final states
      iqx=0
      do iq=iq_start,iq_end
      iqx=iqx+1
      iwx=0
      do iw=iw_start,iw_end
      iwx=iwx+1
      yin(iwx)=-dble(zampa(iq,iw,ip,idx,ipw,igmb,ils))
      yin2(iwx)=-imag(zampa(iq,iw,ip,idx,ipw,igmb,ils))
      end do
      call spline (4, xin, yin, 1,xout,yout)
      call spline (4, xin, yin2, 1,xout,yout2)
      yin3(iqx)=yout(1)
      yin4(iqx)=yout2(1)
      end do
      call spline (4, xin2, yin3, 1,xout2,yout)
      call spline (4, xin2, yin4, 1,xout2,yout2)
      zzz=(yout(1)+zei*yout2(1))*iax
      zmtx(idx,ipw,igmb,ils)=zzz
      zmtx(idxx,ipw,igmb,ils)=zzz*pha
  300 continue 
      endif                     ! tpin.ge.abs(tpinz)
  310 continue
! adding pion pole term
      if(i_pion_pole==1.and.mode.gt.0)then
        fpio2=fpio**2
        qc0=(wcm**2-fnuci**2-Q2)/(2*wcm)
        qc2=Q2+qc0**2
        qc=sqrt(qc2)
        fac=1/(-Q2-fpio2)
      do 550 ipw=1,njLs
      do 550 igmb=igmb_final_meson,igmb_final_meson
      ic        = icgmb(igmb)
      do 550 ils=1,nLsdt(ic,ipw)
      zp=(qc0*zmtx(3,ipw,igmb,ils)
     &    -qc*zmtx(7,ipw,igmb,ils))*fac
      zm=(qc0*zmtx(4,ipw,igmb,ils)
     &    -qc*zmtx(8,ipw,igmb,ils))*fac
      zmtx(3,ipw,igmb,ils)=zmtx(3,ipw,igmb,ils)-qc0*zp
      zmtx(7,ipw,igmb,ils)=zmtx(7,ipw,igmb,ils)-qc*zp
      zmtx(4,ipw,igmb,ils)=zmtx(4,ipw,igmb,ils)-qc0*zm
      zmtx(8,ipw,igmb,ils)=zmtx(8,ipw,igmb,ils)-qc*zm
  550 continue
      endif                     ! i_pion_pole
      else                      ! mode=10 electron scattering; mode=20 photon absorption
      zmtx=0
      endif                     ! mode
! vector current interpolation -> V-A
      if(itiz==1)then         ! proton target
      do 610 ipw=1,njLs
      tpin=dble(itpind(ipw))/2
      if(tpin+eps_tpin.gt.abs(tpinz).and.tmax.gt.tpin)then
        phv=(-1)**((jpind(ipw)-1)/2 + Lpind(ipw)/2+1)
        if(jpind(ipw)==1)then
          idxp_start=2
        else
          idxp_start=1  
        endif
      do 600 igmb=igmb_final_meson,igmb_final_meson
      ic        = icgmb(igmb)
      do 600 ils=1,nLsdt(ic,ipw)
      do idxp=idxp_start,idxp_mx_v
        idx=id1(idxp)
        idxx=id2(idxp)
! stable final states        
      iqx=0
      do iq=iq_start,iq_end
      iqx=iqx+1
      iwx=0
      do iw=iw_start,iw_end
      iwx=iwx+1
      yin(iwx)=dble(zampv(iq,iw,ip,idx,ipw,igmb,ils))
      yin2(iwx)=imag(zampv(iq,iw,ip,idx,ipw,igmb,ils))
      end do
      call spline (4, xin, yin, 1,xout,yout)
      call spline (4, xin, yin2, 1,xout,yout2)
      yin3(iqx)=yout(1)
      yin4(iqx)=yout2(1)
      end do
      call spline (4, xin2, yin3, 1,xout2,yout)
      call spline (4, xin2, yin4, 1,xout2,yout2)
      zzz=(vfac*(yout(1)+zei*yout2(1)))*ivec
      zmtx(idx,ipw,igmb,ils)=zmtx(idx,ipw,igmb,ils)
     &  +zzz
      zmtx(idxx,ipw,igmb,ils)=zmtx(idxx,ipw,igmb,ils)
     &  +zzz*phv
      if(idxp==3)then
! zth component from 0-th component using current conservation
      zmtx(7,ipw,igmb,ils)=zmtx(7,ipw,igmb,ils)+zzz*xxx
      zmtx(8,ipw,igmb,ils)=zmtx(8,ipw,igmb,ils)+zzz*xxx*phv
      endif 
      end do                    ! idxp
  600 continue 
      endif                     ! tpin.ge.abs(tpinz)
  610 continue

      else if(itiz==-1)then      ! neutron target
      do 615 ipw=1,njLs
        phv=(-1)**((jpind(ipw)-1)/2 + Lpind(ipw)/2+1)
      tpin=dble(itpind(ipw))/2
        if(jpind(ipw)==1)then
          idxp_start=2
        else
          idxp_start=1  
        endif 
      if((mode.lt.10.or.itpind(ipw)==3).and.tmax.gt.tpin
     &    .and.tpin+eps_tpin.gt.abs(tpinz))then ! weak or (EM and isospin=3/2)
      do 605 igmb=igmb_final_meson,igmb_final_meson
      ic        = icgmb(igmb)
      do 605 ils=1,nLsdt(ic,ipw)
      do idxp=idxp_start,idxp_mx_v
        idx=id1(idxp)
        idxx=id2(idxp)
! stable final states        
      iqx=0
      do iq=iq_start,iq_end
      iqx=iqx+1
      iwx=0
      do iw=iw_start,iw_end
      iwx=iwx+1
      yin(iwx)=dble(zampv(iq,iw,ip,idx,ipw,igmb,ils))
      yin2(iwx)=imag(zampv(iq,iw,ip,idx,ipw,igmb,ils))
      end do
      call spline (4, xin, yin, 1,xout,yout)
      call spline (4, xin, yin2, 1,xout,yout2)
      yin3(iqx)=yout(1)
      yin4(iqx)=yout2(1)
      end do
      call spline (4, xin2, yin3, 1,xout2,yout)
      call spline (4, xin2, yin4, 1,xout2,yout2)
      zzz=(vfac*(yout(1)+zei*yout2(1)))*ivec
      zmtx(idx,ipw,igmb,ils)=zmtx(idx,ipw,igmb,ils)
     &  +zzz
      zmtx(idxx,ipw,igmb,ils)=zmtx(idxx,ipw,igmb,ils)
     &  +zzz*phv
      if(idxp==3)then
! zth component from 0-th component using current conservation
      zmtx(7,ipw,igmb,ils)=zmtx(7,ipw,igmb,ils)+zzz*xxx
      zmtx(8,ipw,igmb,ils)=zmtx(8,ipw,igmb,ils)+zzz*xxx*phv
      endif 
      end do                    ! idxp
  605 continue 
      else if(itpind(ipw)==1.and.tpin+eps_tpin.gt.abs(tpinz))then ! EM and I=1/2
      do 606 igmb=igmb_final_meson,igmb_final_meson
      ic        = icgmb(igmb)
      do 606 ils=1,nLsdt(ic,ipw)
      do idxp=idxp_start,idxp_mx_v
        idx=id1(idxp)
        idxx=id2(idxp)
! stable final states        
      iqx=0
      do iq=iq_start,iq_end
      iqx=iqx+1
      iwx=0
      do iw=iw_start,iw_end
      iwx=iwx+1
      yin(iwx)=dble(zampv_is(iq,iw,ip,idx,ipw,igmb,ils))
      yin2(iwx)=imag(zampv_is(iq,iw,ip,idx,ipw,igmb,ils))
      end do
      call spline (4, xin, yin, 1,xout,yout)
      call spline (4, xin, yin2, 1,xout,yout2)
      yin3(iqx)=yout(1)
      yin4(iqx)=yout2(1)
      end do
      call spline (4, xin2, yin3, 1,xout2,yout)
      call spline (4, xin2, yin4, 1,xout2,yout2)
      zzz=(vfac*(yout(1)+zei*yout2(1)))*ivec
      zmtx(idx,ipw,igmb,ils)=zmtx(idx,ipw,igmb,ils)
     &  +zzz
      zmtx(idxx,ipw,igmb,ils)=zmtx(idxx,ipw,igmb,ils)
     &  +zzz*phv
      if(idxp==3)then
! zth component from 0-th component using current conservation
      zmtx(7,ipw,igmb,ils)=zmtx(7,ipw,igmb,ils)+zzz*xxx
      zmtx(8,ipw,igmb,ils)=zmtx(8,ipw,igmb,ils)+zzz*xxx*phv
      endif 
      end do                    ! idxp
  606 continue 
      endif 
  615 continue
      endif                     ! itiz

! isoscalar vector current for N!      

      if(mode.le.-1)then
      do 720 ipw=1,njLs
      tpin=dble(itpind(ipw))/2
      if(itpind(ipw)==1.and.tpin+eps_tpin.gt.abs(tpinz))then ! I=1/2
        phv=(-1)**((jpind(ipw)-1)/2 + Lpind(ipw)/2+1)
        if(jpind(ipw)==1)then
          idxp_start=2
        else
          idxp_start=1  
        endif 
      do 700 igmb=igmb_final_meson,igmb_final_meson
      ic        = icgmb(igmb)
      do 700 ils=1,nLsdt(ic,ipw)
      do idxp=idxp_start,idxp_mx_v
        idx=id1(idxp)
        idxx=id2(idxp)
      iqx=0
      do iq=iq_start,iq_end
      iqx=iqx+1
      iwx=0
      do iw=iw_start,iw_end
      iwx=iwx+1
      yin(iwx)=dble(zampv_is(iq,iw,ip,idx,ipw,igmb,ils))
      yin2(iwx)=imag(zampv_is(iq,iw,ip,idx,ipw,igmb,ils))
      end do
      call spline (4, xin, yin, 1,xout,yout)
      call spline (4, xin, yin2, 1,xout,yout2)
      yin3(iqx)=yout(1)
      yin4(iqx)=yout2(1)
      end do
      call spline (4, xin2, yin3, 1,xout2,yout)
      call spline (4, xin2, yin4, 1,xout2,yout2)
      zzz=vvfac(itiz)*(yout(1)+zei*yout2(1))*ivec
      zmtx(idx,ipw,igmb,ils)=zmtx(idx,ipw,igmb,ils)
     &  +zzz
      zmtx(idxx,ipw,igmb,ils)=zmtx(idxx,ipw,igmb,ils)
     &  +zzz*phv
      if(idxp==3)then
! zth component from 0-th component using current conservation
      zmtx(7,ipw,igmb,ils)=zmtx(7,ipw,igmb,ils)+zzz*xxx
      zmtx(8,ipw,igmb,ils)=zmtx(8,ipw,igmb,ils)+zzz*xxx*phv
      endif 
      end do                    ! idxp
  700 continue 
      endif 
  720 continue
      endif                     ! mode.le.-1 
      return

      end




      
      subroutine ylmsub(lmax,z,bleg)
      implicit real*8(a-h,o-z)
      common / cscal / pi,fm
      dimension bleg(0:8,-8:8),bc(0:30),bb(0:30,0:30)
      save bb,iniylm

      if(iniylm.eq.0)then
      bc(0)   = 1
      do 100 k= 1,30
      bc(k)   = bc(k-1)*dble(k)
  100 continue
      do 200 k1= 0,30
      do 200 k2= 0,30
      bb(k1,k2)= bc(k1)/bc(k2)
  200 continue
      iniylm   = 1
      end if

      z1         = sqrt(1.d0 - z**2) + 1d-20 ! avoid 0; corrected on 030916
!      z1         = sqrt(1.d0 - z**2)
      z2         = z/z1
      bleg(0,0)  = 1.d0
      do 300 l   = 1,lmax
      bleg(l,l)  = z1**l/dble(2**l)*bb(2*l,l)
      bleg(l,l-1)= z2*bleg(l,l)
      if(l.eq.1) go to 300
      do 310 m   = l-2,0,-1
      bleg(l,m)  = (-bleg(l,m+2) + dble(2*(m+1))*z2*bleg(l,m+1))
     &              / dble((l-m)*(l+m+1))
  310 continue
  300 continue

      do 400 l   = 0,lmax
      do 400 m   = 0,l
      fac        = sqrt(dble(2*l+1)/4.d0/pi*bb(l-m,l+m))*
     &             dble((-1)**m)
      bleg(l,m)  = bleg(l,m)*fac
      bleg(l,-m) = bleg(l,m)*dble((-1)**m)
  400 continue
      return
      end
      

      subroutine lorentz_trans(i,j,pcm,xlr)
      implicit real*8 (a-h,o-y)
      implicit complex*16 (z)
      dimension pcm(0:3),xlr(0:3,0:3),qlab(0:3),q2cm(0:3)
     &  ,xlrx(0:3,0:3),rot(3,3)
      save rot
      xm=sqrt(pcm(0)**2-pcm(1)**2-pcm(2)**2-pcm(3)**2)
      ifac=(-1)**i              ! ifac=-1 for "Lab to CM"; =1 for "CM to Lab"
      rxm=1/xm
      xlr(0,1)=pcm(1)*rxm*ifac
      xlr(0,2)=pcm(2)*rxm*ifac
      xlr(0,3)=pcm(3)*rxm*ifac
      xlr(0,0)=pcm(0)*rxm
      xxx=1/(xm*(xm+pcm(0)))
      xlr(1,1)=1+pcm(1)**2*xxx
      xlr(1,2)=pcm(1)*pcm(2)*xxx
      xlr(1,3)=pcm(1)*pcm(3)*xxx
      xlr(2,2)=1+pcm(2)**2*xxx
      xlr(2,3)=pcm(2)*pcm(3)*xxx
      xlr(3,3)=1+pcm(3)**2*xxx
      do n=0,3
      do m=n+1,3
        xlr(m,n)=xlr(n,m)
      end do
      end do
      
      return
      end
!------------------------------------------------------------
      subroutine bifc
!------------------------------------------------------------
      implicit real*8(a-h,o-z)
      parameter (n=100,m=50)
      common / fdbn / h(0:n),dh(-1:m),bb(0:n,0:n)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      e=1.0d0 
      h(0)=e 
      dh(-1)=e 
      bb(0,0)=e 
      sd=e 
      do 10 i=1,n 
        s=i 
        sd=sd*sqrt(s) 
   10 h(i)=sd 
      sd=e 
      do 20 i=0,m 
        s=i+i+1 
        sd=sd*sqrt(s) 
   20 dh(i)=sd
      do 30 i=1,n 
        bb(i,0)=e 
        bb(0,i)=e 
        bb(i,i)=e 
        s=e
        do 31 j=1,i/2 
          l=i-j 
          s=s*dble(l+1)/dble(j)
          sd=sqrt(s) 
          bb(j,i)=sd 
          bb(l,i)=sd 
          bb(i,j)=s 
          bb(i,l)=s
   31   continue
   30 continue
      return 
      end
!------------------------------------------------------------
      double precision function cbg(a,x,b,y,c,z)
!------------------------------------------------------------
      implicit real*8(a-h,o-z)
      common / fdbn / h(0:100),dh(-1:50),bb(0:100,0:100)
      parameter (de=0.01d0,t=1.0d0) 
      prt(i)=1-2*mod(i,2)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      cbg=0.d0 
      if(max(abs(a-b)-c,c-a-b,abs(x+y-z),
     &  x-a,-x-a,y-b,-y-b).gt.de) return  
      a1=a+de 
      b1=b+de 
      c1=c+de 
      ja=int(a1-x)
      jb=int(b1+y)
      jc=int(c1+z)
      ka=int(a1+a)
      kb=int(b1+b)
      kc=int(c1+c)
!      ks=int(ka+kb)
      is=int(a1+b+c)
      ia=is-ka 
      ib=is-kb 
      ic=is-kc
      lmin=max(0,ja-ib,jb-ia) 
      lmax=min(ic,ja,jb)
      fl=-prt(lmin) 
      do 10 l=lmin,lmax 
        fl=-fl
   10 cbg=cbg+fl*bb(ic,l)*bb(ib,ja-l)*bb(ia,jb-l)
      cbg=cbg*sqrt((dble(kc)+t)/(dble(is)+t))*bb(ib,kc)*bb(ic,ka)
     &  /(bb(kb,is)*bb(ja,ka)*bb(jb,kb)*bb(jc,kc))
      return 
      end
!-----------------------------------------------------------
      subroutine spline (n, x, y, nout,xout,yout)
      implicit real*8(a-h,o-z)
      parameter(ndim=60)
      dimension x(ndim),y(ndim),b(ndim),c(ndim),d(ndim),
     &          xout(ndim),yout(ndim)

!  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!  for a cubi! interpolating spline

!    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3

!    for  x(i) .le. x .le. x(i+1)

!  input..

!    n = the number of data points or knots (n.ge.2)
!    x = the abscissas of the knots in strictly increasing order
!    y = the ordinates of the knots

!  output..

!    b, c, d  = arrays of spline coefficients as defined above.

!  using  p  to denote differentiation,

!    y(i) = s(x(i))
!    b(i) = sp(x(i))
!    c(i) = spp(x(i))/2
!    d(i) = sppp(x(i))/6  (derivative from the right)

!  the accompanying function subprogram  seval  can be used
!  to evaluate the spline.


!      integer nm1, ib, i
!      double precision t

      nm1 = n-1
      if ( n .lt. 2 ) return
      if ( n .lt. 3 ) go to 50

!  set up tridiagonal system

!  b = diagonal, d = offdiagonal, ! = right hand side.

      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do 10 i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
   10 continue

!  end conditions.  third derivatives at  x(1)  and  x(n)
!  obtained from divided differences

      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.d0
      c(n) = 0.d0
      if ( n .eq. 3 ) go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))

!  forward elimination

   15 do 20 i = 2, n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
   20 continue

!  back substitution

      c(n) = c(n)/b(n)
      do 30 ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
   30 continue

!  c(i) is now the sigma(i) of the text

!  compute polynomial coefficients

      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2*c(n))
      do 40 i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3*c(i)
   40 continue
      c(n) = 3*c(n)
      d(n) = d(n-1)
      go to 100

   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.d0
      d(1) = 0.d0
      b(2) = b(1)
      c(2) = 0.d0
      d(2) = 0.d0
      go to 100

!    calculate yout

 100  continue
      do 200 j=1,nout
      xo      = xout(j)
      i       = 1
 210  continue
      if(x(i).le.xo.and.x(i+1).ge.xo) go to 220
      i = i+1
      if(i.gt.n) go to 230
      go to 210
  220 del = xo - x(i)
      yout(j) = y(i)+ b(i)*del + c(i)*del**2+d(i)*del**3
      go to 200
  230 write(*,*)' spline interpolation out of range!!'
  200 continue
      return
      end

      subroutine rspin(ptot,p00,pc,zmx)
      implicit real*8(a-h,o-y)
      implicit complex*16(z)
      common / ccoup / fpio,fnuc,flep,fnuci,fdeu,feta,fmfin
      common / cscal / pi,fm
      dimension p00(0:3),pc(0:3), ptot(0:3),pcx(0:3)
      dimension zmx(2,2),zmxs(2,2),
     1 zmx1(2,2),zmxp(2,2),zmxpc(2,2),zmxptot(2,2)
      dimension zsigm(2,2,0:3)
      dimension pv(3),pcv(3),vbeta(3)
      data idata / 1 /
      save idata,zsigm

!      Pauli matrices

      amn=fnuc
       if(idata==1)then
       Z=cmplx(0.,1.)
       DO 11 I1=1,2
       DO 11 I2=1,2
       DO 11 IA=0,3
   11  ZSIGM(I1,I2,IA)=0.
       ZSIGM(1,1,0)=1
       ZSIGM(2,2,0)=1
       ZSIGM(1,2,1)=1
       ZSIGM(2,1,1)=1
       ZSIGM(1,2,2)=-Z
       ZSIGM(2,1,2)=Z
       ZSIGM(1,1,3)=1
       ZSIGM(2,2,3)=-1
       idata=0 
       endif 

       do i=1,3
       vbeta(i)=ptot(i)/ptot(0)
       end do

        p0=p00(0)
        amtot=sqrt(ptot(0)**2-ptot(1)**2-ptot(2)**2-ptot(3)**2)
        do  i1=1,2
        do  i2=1,2
        zmxp(i1,i2)=(p0+amn)*zsigm(i1,i2,0)
        zmxpc(i1,i2)=(pc(0)+amn)*zsigm(i1,i2,0)
        zmxptot(i1,i2)=(ptot(0)+amtot)*zsigm(i1,i2,0)
        do i=1,3
        zmxp(i1,i2)=zmxp(i1,i2)+p00(i)*zsigm(i1,i2,i)
        zmxpc(i1,i2)=zmxpc(i1,i2)-pc(i)*zsigm(i1,i2,i)
        zmxptot(i1,i2)=zmxptot(i1,i2)-ptot(i)*zsigm(i1,i2,i)
        end do
        end do
        end do
        fac=1d0/sqrt(2d0*amn*(p0+amn))
        zmxp=zmxp*fac
        fac=1d0/sqrt(2d0*amn*(pc(0)+amn))
        zmxpc=zmxpc*fac
        fac=1d0/sqrt(2d0*amtot*(ptot(0)+amtot))
        zmxptot=zmxptot*fac

        do i1=1,2
        do i2=1,2
        zmx1(i1,i2)=0.
        do i3=1,2
        zmx1(i1,i2)=zmx1(i1,i2)+zmxptot(i1,i3)*zmxp(i3,i2)
        end do
        end do
        end do

        do i1=1,2
        do i2=1,2
        zmx(i1,i2)=0.
        do i3=1,2
        zmx(i1,i2)=zmx(i1,i2)+zmxpc(i1,i3)*zmx1(i3,i2)
        end do
        end do
        end do
!        do i1=1,2
!        do i2=1,2
!        zmxs(i1,i2)=conjg(zmx(i1,i2))
!        end do
!        end do
       return
       end
      subroutine setdfun(x,dfun)
      implicit real*8(a-h,o-z)
      parameter (njmx=5)
      common /para_dfn / jmax
      dimension dfun(2*njmx-1,-2*njmx+1:2*njmx-1,-2*njmx+1:2*njmx-1)
!      parameter(max_l=27)
!      common / cdff / xgau(100),wgau(100),dfun(2*njmx-1,-5:5,-5:5,100)
!     & ,fleg(0:max_l,100)
!      common / cdfi / meshx,mxx,mxj,mxm
!      dimension fle(0:max_l),fled(-1:max_l),fledd(-1:max_l)
!       dimension xg(3),wg(3)
!       data xg/ -0.7745966692D0,0.d0,0.7745966692D0/
!       data wg/ 0.5555555555d0,0.8888888888d0,0.5555555555d0/
!       data ngaus/3/

      dfun = 0
!      fleg = 0

!      mxx       = meshx * ngaus
!      dgux      = 2.d0/dble(2*meshx)
!      idx       = 1
!      do 110 nx = 1,meshx
!      do 110 ng = 1,ngaus
!      xgau(idx)   = dgux*(xg(ng)+dble(2*nx-1))-1.d0
!      wgau(idx)   = wg(ng)*dgux
!      x         = xgau(idx)

!       call legen(x,fle)
!       do 111 lx = 0,max_l
!  111  fleg(lx,idx) = fle(lx)

      ss        = sqrt((1.d0 - x)/2.d0)
      cc        = sqrt((1.d0 + x)/2.d0)

      do 120 lx = 1,jmax,2
!      do 120 lx = 1,2*njmx-1,2
!      maxm      = min(lx,mxm)
      maxm      = lx
      do 130 mf = -maxm,maxm,2
      do 140 mi = -maxm,maxm,2
      jj        = (lx - 1)/2
      mfm       = (mf - 1)/2
      mfp       = (mf + 1)/2
      mim       = (mi - 1)/2
      mip       = (mi + 1)/2
      df1       = 0
      df2       = 0
      df3       = 0
      df4       = 0
      if(jj.ge.abs(mfm).and.jj.ge.abs(mim))then
      fac       = sqrt(dble((lx + mf)*(lx + mi)))/dble(2*lx)
      df1       = fblmmx(jj,mfm,mim,cc,ss)*fac*cc
      end if

      if(jj.ge.abs(mfp).and.jj.ge.abs(mip))then
      fac       = sqrt(dble((lx - mf)*(lx - mi)))/dble(2*lx)
      df2       = fblmmx(jj,mfp,mip,cc,ss)*fac*cc
      end if

      if(jj.ge.abs(mfm).and.jj.ge.abs(mip))then
      fac       = sqrt(dble((lx + mf)*(lx - mi)))/dble(2*lx)
      df3       =-fblmmx(jj,mfm,mip,cc,ss)*fac*ss
      end if

      if(jj.ge.abs(mfp).and.jj.ge.abs(mim))then
      fac       = sqrt(dble((lx - mf)*(lx + mi)))/dble(2*lx)
      df4       = fblmmx(jj,mfp,mim,cc,ss)*fac*ss
      end if

      dfun(lx,mf,mi) = df1 + df2 + df3 + df4

 140  continue
 130  continue
 120  continue

!      idx       = idx + 1
! 110  continue

      return
      end
      real*8 function fblmmx(l,mf,mi,cc,ss)
      implicit real*8(a-h,o-z)
      parameter (n=100,m=50)
      common / fdbn / h(0:n),dh(-1:m),bb(0:n,0:n)

      fblmmx  = 0
      if(l.lt.0.or.abs(mf).gt.l.or.abs(mi).gt.l)return

!      t2     = theta/2.d0
!      c!     = cos(t2)
!      ss     = sin(t2)
      jmip   = l  + mi
      jmim   = l  - mi
      jmfp   = l  + mf
      jmfm   = l  - mf
      mfmim  = mf - mi
      iicos  = 2*l - mfmim
      iisin  = mfmim
      kmax   = min(jmip,jmfm)
      kmin   = max(0,-mfmim)
      if(kmax.lt.kmin)return
      sum    = 0
      do 100 kx = kmin,kmax
      phase     = (-1)**(mfmim+kx)
      factor    = h(jmip)*h(jmim)*h(jmfp)*h(jmfm)
     &           /(h(jmip-kx)*h(kx)*h(jmfm-kx)*h(kx+mfmim))**2
      sum       = phase*factor*cc**(iicos-2*kx)*ss**(iisin+2*kx)
     &           + sum
 100  continue
      fblmmx     = sum
      return
      end function

