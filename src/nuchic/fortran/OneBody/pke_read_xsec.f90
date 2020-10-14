 module quasi_el
    use libutilities
    use libinterpolate
    implicit none
    real*8 :: mp,mn,mqe
    real*8 :: hbarc
    real*8, parameter :: pi= 4.0d0*atan(1.0d0),c0=1.0d0/16.0d0/pi/pi
    real*8, allocatable, save ::pke_p(:,:),xe_p(:),dp_p(:)
    real*8, allocatable, save ::pke_n(:,:),xe_n(:),dp_n(:)
    real*8, allocatable, save :: p(:)
    
    integer*4, save, private :: fg,nZ,nA,iform
    integer*4, save :: np,nen,nep
    real*8, save, private :: kF

    type(interp2d) :: pke_p_interp, pke_n_interp

    contains


    subroutine init_pke(fname_pkep,fname_pken,fg_in,nZ_in,nA_in,kF_in,iform_in)
        implicit none
        integer*4 :: fg_in,i,j,nZ_in,nA_in,iform_in
        real*8 :: he_p,he_n
        real*8 :: norm,pmax, hp,kF_in
        character*50 :: fname_pkep,fname_pken
        call init(constants)
        mp = constants%mp
        mn = constants%mn
        mqe = constants%mqe
        hbarc = constants%hbarc
        fg=fg_in
        nZ=nZ_in
        nA=nA_in
        kF=kF_in
        iform= iform_in
        if(fg.ne.1)then
          open(unit=4,file=fname_pkep,status='unknown',form='formatted')
          read(4,*) nep,np
          allocate(p(np),pke_p(nep,np),dp_p(np),xe_p(nep))
          do j=1,np
             read(4,*) p(j)
             read(4,'(4(f6.1,2x,e10.3))')(xe_p(i),pke_p(i,j),i=1,nep)
          enddo
          close(4)
          pke_p=pke_p*(xe_p(2)-xe_p(1))
          pke_p_interp = interp2d(xe_p, p, pke_p, nep, np)
!
          open(unit=4,file=fname_pken,status='unknown',form='formatted')
          read(4,*) nen,np
          allocate(pke_n(nen,np),dp_n(np),xe_n(nen))
          do j=1,np
             read(4,*) p(j)
             read(4,'(4(f6.1,2x,e10.3))')(xe_n(i),pke_n(i,j),i=1,nen)
          enddo
          close(4)
          pke_n=pke_n*(xe_n(2)-xe_n(1))
          pke_n_interp = interp2d(xe_n, p, pke_n, nep, np)

          !... change dimension to p[MeV] and [MeV**-4]
          !p=p*hbarc
          !Pke_p=Pke_p/hbarc**3
          !Pke_n=Pke_n/hbarc**3       
          pmax=p(np)
          hp=p(2)-p(1)!pmax/dble(nbox)
       else
          write(6,*) 'we did not code the FG case for the Asymmetric nuclei'
       endif
    
       norm=0.0d0
       do j=1,np
          dp_p(j)=sum(pke_p(:,j))!*he_p
       enddo
    
       norm=sum(p(:)**2*dp_p(:))*4.0d0*pi*hp
       pke_p=pke_p/norm
       write(6,*)'n(k) norm initial for protons=', norm
    
       norm=0.0d0
       do j=1,np
          dp_n(j)=sum(pke_n(:,j))!*he_n
       enddo
       norm=sum(p(:)**2*dp_n(:))*4.0d0*pi*hp
       pke_n=pke_n/norm
       write(6,*)'n(k) norm initial for neutrons=', norm

end subroutine



  subroutine f_eval(in,p_4,pf_4,e,mom,w,qval,thetalept,ee,f_o)
    use mathtool
    implicit none
    real*8, parameter :: eps=5.0d0,small=1e-15
    real*8 :: e, mom
    
    integer*4 :: in
    real*8 :: p_4(4),pf_4(4)
    real*8 :: w,wt,qval,thetalept,ee,f_o,pke,xp,xpf
    real*8 :: sig,arg,delta_w

    f_o=0.0d0
    if(in.eq.1) then
        pke = pke_p_interp%call(e, mom)
    elseif(in.eq.2) then
        pke = pke_n_interp%call(e, mom)
    endif       
    xp=sqrt(sum(p_4(2:4)**2))
    xpf=sqrt(sum(pf_4(2:4)**2))
    wt=w-abs(e)+mqe-p_4(1)
    

    arg=wt+p_4(1)-pf_4(1)
    delta_w=fdelta(arg,eps)
    if (delta_w.gt.small)then
         call cc1(in,qval/hbarc,w,wt,xp/hbarc,xpf/hbarc,p_4/hbarc,pf_4/hbarc,ee,thetalept,iform,sig)
         f_o=xp**2*pke*(dble(nZ)*sig)*2.0d0*pi*delta_w*2.0d0
    endif
    
    return
  end subroutine f_eval


  end module quasi_el
