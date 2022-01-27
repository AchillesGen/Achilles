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
    
    integer*4, save, private :: fg,iform
    integer*4, save :: np,nen,nep

    type(interp2d) :: pke_p_interp, pke_n_interp

    contains

    subroutine delete_pke()
        implicit none
        if(fg.ne.1) then
            call delete(pke_p_interp)
            call delete(pke_n_interp)
            deallocate(pke_p,dp_p,xe_p)
            deallocate(p,pke_n,dp_n,xe_n)
        endif
    end subroutine

    subroutine init_pke(fname_pkep,fname_pken,fg_in,iform_in)
        implicit none
        integer*4 :: fg_in,i,j,iform_in
        real*8 :: he_p,he_n
        real*8 :: norm,pmax, hp
        real*8 :: dpe, ee, pp, pke
        character(len=*) :: fname_pkep,fname_pken
        call init(constants)
        mp = constants%mp
        mn = constants%mn
        mqe = constants%mqe
        hbarc = constants%hbarc
        iform= iform_in
        fg = fg_in
        if(fg.ne.1)then
          open(unit=4,file=fname_pkep,status='unknown',form='formatted', action='read')
          read(4,*) nep,np
          allocate(p(np),pke_p(nep,np),dp_p(np),xe_p(nep))
          do j=1,np
             read(4,*) p(j)
             read(4,'(4(f6.1,2x,e10.3))')(xe_p(i),pke_p(i,j),i=1,nep)
          enddo
          close(4)
          norm=0.0d0
          hp=p(2)-p(1)
          he_p=xe_p(2)-xe_p(1)
          do j=1,np
             dp_p(j)=sum(pke_p(:,j))*he_p
          enddo
    
          norm=sum(p(:)**2*dp_p(:))*4.0d0*pi*hp
          pke_p_interp = interp2d(xe_p, p, pke_p, nep, np, 1)
          write(6,*)'n(k) norm initial for protons=', norm

          open(unit=4,file=fname_pken,status='unknown',form='formatted', action='read')
          read(4,*) nen,np
          allocate(pke_n(nen,np),dp_n(np),xe_n(nen))
          do j=1,np
             read(4,*) p(j)
             read(4,'(4(f6.1,2x,e10.3))')(xe_n(i),pke_n(i,j),i=1,nen)
          enddo
          close(4)
          norm=0.0d0
          he_n=xe_n(2)-xe_n(1)

          do j=1,np
             dp_n(j)=sum(pke_n(:,j))*he_n
          enddo
          norm=sum(p(:)**2*dp_n(:))*4.0d0*pi*hp
          pke_n_interp = interp2d(xe_n, p, pke_n, nep, np, 1)
          write(6,*)'n(k) norm initial for neutrons=', norm

     !     do i=1,nep
     !         dpe=0.0d0
     !         do j=1,np
     !            dpe=dpe+pke_p(i,j)*p(j)**2*hp
     !          enddo   
     !          write(1001,*) xe_p(i), dpe
     !     enddo  



    !      he_p=(xe_p(nep)-xe_p(1))/(2.*nep)
    !      hp=(p(np)-p(1))/(2.*np)
    !      do i=1,2*nep
    !          ee=(dble(i)-0.5d0)*he_p
    !          dpe=0.0d0
    !          do j=1,np*2
    !             pp=(dble(j)-0.5d0)*hp
    !             pke = pke_p_interp%call(ee, pp)
    !             dpe=dpe+pke*pp**2*hp
    !           enddo   
    !           write(1000,*) ee, dpe
    !      enddo     
    !   stop
       endif
    end subroutine

    subroutine f_eval(p_4,pf_4,k_4,kp_4,e,mom,w,qval,thetalept,ee,nZ,nA,kF,fp_o, fn_o)
        use mathtool
        implicit none
        real*8, parameter :: eps=5.0d0,small=1e-15
        real*8, intent(in) :: e, mom
        integer, intent(in) :: nZ, nA 
        real*8, intent(in) :: kF
        real*8 :: p_4(4)
        real*8, intent(in) :: pf_4(4),k_4(4),kp_4(4)
        real*8, intent(in) :: w, qval, thetalept, ee
        real*8 :: wt, pkep, pken, xp, xpf, pke
        real*8, intent(out) :: fp_o, fn_o
        real*8 :: sig_p, sig_n, arg, delta_w, norm

        fp_o=0.0d0      
        fn_o=0.0d0      
        xp=sqrt(sum(p_4(2:4)**2))
        xpf=sqrt(sum(pf_4(2:4)**2))

        if(fg.ne.1) then
           pkep = pke_p_interp%call(e, mom)
           pken = pke_n_interp%call(e, mom)
           p_4(1)=sqrt(xp**2+mqe**2)
           wt=w-abs(e)+mqe-p_4(1)
           if (wt.lt.0.0d0) return
           call cc1(qval/hbarc,w,wt,xp/hbarc,xpf/hbarc,p_4/hbarc,pf_4/hbarc,k_4/hbarc,kp_4/hbarc,thetalept,iform,sig_p,sig_n)
           ! f_o=xp**2*pke*(dble(nZ)*sig)*2.0d0*pi*delta_w*2.0d0
           fp_o=pkep*sig_p/dble(nZ)
           ! fp_o=sig_p
           fn_o=pken*sig_n/dble(nA-nZ)
        else
            if(xp.gt.kF) then
                fp_o=0.0d0   
                fn_o=0.0d0
                return
            endif
            norm=4.0d0*pi/3.0d0*kF**3  
            pke=1.0/norm   
            p_4(1)=sqrt(xp**2+mqe**2)
            wt=w-abs(e)
            call cc1(qval/hbarc,w,wt,xp/hbarc,xpf/hbarc,p_4/hbarc,pf_4/hbarc,k_4/hbarc,kp_4/hbarc,thetalept,iform,sig_p, sig_n)
            ! f_o=xp**2*pke*(dble(nZ)*sig)*2.0d0*pi*delta_w*2.0d0
            fp_o=pke*sig_p*pf_4(1)/xp/qval
            fn_o=pke*sig_n*pf_4(1)/xp/qval
        endif
    
        return
    end subroutine f_eval

end module quasi_el
