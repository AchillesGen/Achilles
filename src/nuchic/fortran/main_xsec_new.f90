program response_ia
  use mathtool
  use mympi
  use dirac_matrices
  implicit none
    integer*4, parameter :: nev=1000000,neq=100,nvoid=10
    integer*8, allocatable :: irn(:),irn0(:)
    real*8, parameter :: mp=938.272046d0,mn=939.56563d0,hbarc=197.327053d0
    real*8, parameter :: xme=0.511d0
    real*8, parameter :: small=1e-12


    integer*4 :: np0,nh0,nA,nZ,np,ne,nc,nw,i,j,k,nskip,iform,ie,iemax,nbox,fg
    integer*4, allocatable :: ip_o(:),ip_n(:),ie_o(:),ie_n(:)
    integer*4 :: nwlk
    integer*4 :: i_acc,i_avg,i_acc_tot,i_avg_tot,iv,iw
    real*8, allocatable :: p0(:),pke0(:,:),p(:),dp(:),cost_d(:)
    real*8, allocatable :: w(:),wt(:),sigma(:),Q2(:),pke(:,:),xe(:),he
    real*8 :: ee,thetalept,coste,eef,sig_p,sig_n, cost_te,norm,mqef,dummy
    real*8 :: pmax,hp,wmax,hw,pi,rdummy,mspec,mnuc,mqe,qval,qfm,eps,sig,arg,mstar
    real*8 :: cost,cost_rel,pf,dpf,kf,TA,dwfold,ffold,hwfold,hc,delta_w,espec,epf,ep
    real*8 :: r_avg,r_err,r_avg_tot,r_err_tot,mpi
    real*8, allocatable :: g_o(:),g_n(:),f_o(:),f_n(:)
    real*8 :: ti,tf
    character * 40 :: nk_fname
    character*50 :: fname,en_char,theta_char,int_char,dec_char
    logical :: explicit_d
    call init0()
! inputs
    if (myrank().eq.0) then
       read(5,*) nwlk
       read(5,*) ee,thetalept
       read(5,*) nZ, nA
       read(5,*) fg
       read(5,*) kF
       read(5,*) nw
       read(5,*) wmax
    endif
    call bcast(nwlk)
    call bcast(ee)
    call bcast(thetalept)
    call bcast(nZ)
    call bcast(nA)
    call bcast(fg)
    call bcast(kF)
    call bcast(nw)
    call bcast(wmax)

    iform=2

!    ti=MPI_Wtime()
    allocate(irn0(nwlk))
    do i=1,nwlk
       irn0(i)=19+i
    enddo
    if (myrank().eq.0) then
       write (6,'(''number of cpus ='',t50,i10)') nproc()
       if (mod(nwlk,nproc()).ne.0) then
          write(6,*)'Error: nwalk must me a multiple of nproc'
          stop
       endif
    endif
    nwlk=nwlk/nproc()

    allocate(irn(nwlk))
    allocate(ip_o(nwlk),ip_n(nwlk),ie_o(nwlk),ie_n(nwlk),g_o(nwlk),g_n(nwlk),f_o(nwlk),f_n(nwlk))
    irn(:)=irn0(myrank()*nwlk+1:myrank()*nwlk+nwlk)

    
    write(int_char,'(i3)') int(thetalept)
    write(dec_char,'(i1)') int(mod(thetalept*10.0d0,10.0d0))

    theta_char=trim(int_char)//'p'//trim(dec_char)
    theta_char=adjustl(theta_char)
    theta_char=trim(theta_char)
    write(en_char,'(i4)') int(ee)
    en_char=adjustl(en_char)
    en_char=trim(en_char)

    fname='C12_'//trim(en_char)//'_'//trim(theta_char)//'.out'
    fname=trim(fname)

! write normalized responses on file
    open(unit=14,file=fname)


    iform=2

! initialize useful constants
    mqe=0.5d0*(mp+mn)
    mqef=mqe/hbarc
    pi=acos(-1.0d0)
    qfm=qval/hbarc
    thetalept=thetalept*pi/180.0d0
    coste= cos(thetalept)
    call dirac_matrices_in(mqef)


! read and interpolate the momentum distribution
    !    open(unit=4,file='sf_O16_e13_hw20_xsec.out',status='unknown',form='formatted')
    if(fg.ne.1)then
       open(unit=4,file='pke/pke12_tot.data',status='unknown',form='formatted')
       read(4,*) ne, np
       allocate(p(np),pke(ne,np),dp(np),xe(ne))
       do j=1,np
          read(4,*) p(j)
          read(4,'(4(f6.1,2x,e10.3))')(xe(i),pke(i,j),i=1,ne)
       enddo
       close(4)
       pmax=p(np)
       hp=p(2)-p(1)!pmax/dble(nbox)
       he=xe(2)-xe(1)
    else
       ne=1
       np=100
       allocate(p(np),pke(ne,np),dp(np),xe(ne))
       p(np)=kF
       hp=p(np)/dble(np)
       xe=0.0d0
       he=1
       do i=1,np
          p(i)=dble(i)*hp
          pke(1,i)=1.0d0
       enddo

    endif
    norm=0.0d0
    do j=1,np
       dp(j)=sum(pke(:,j))*he
    enddo
    norm=sum(p(:)**2*dp(:))*4.0d0*pi*hp
    pke=pke/norm!dble(nZ)
     write(6,*)'n(k) norm initial=', norm
    

! construct omega grid and the form factors
    allocate(w(nw),wt(nw),Q2(nw))
    hw=wmax/dble(nw)
    do i=1,nw
       w(i)=(dble(i)-0.5d0)*hw
       eef = ee - w(i)
       Q2(i) = 2.0d0*ee*eef*(1.0d0 - coste)
    enddo

! compute the response functions in the impulse approximation
    do iw=1,nw
       qval=sqrt( Q2(iw)+ w(iw)**2 )
!
       r_avg=0.0d0
       r_err=0.0d0
       i_acc=0
       i_avg=0
       g_o=0.0d0

       do i=1,nwlk
          call setrn(irn(i))
          do while(g_o(i).le.0.0d0)
             ip_o(i)=1+int(np*ran())
             ie_o(i)=1+int(ne*ran())
             call g_eval(p(ip_o(i)),PkE(ie_o(i),ip_o(i)),g_o(i))
          enddo
          call getrn(irn(i))
       enddo
       
       do iv=1,nev
          do j=1,nwlk
             !call setrn(irn(j))
             ip_n(j)=nint(ip_o(j)+0.05d0*np*(-1.0d0+2.0d0*ran()))
             ie_n(j)=nint(ie_o(j)+0.05d0*ne*(-1.0d0+2.0d0*ran()))
             if (ip_n(j).le.np.and.ip_n(j).ge.1.and.ie_n(j).le.ne.and.ie_n(j).ge.1) then
                call g_eval(p(ip_n(j)),PkE(ie_n(j),ip_n(j)),g_n(j))
             else
                g_n(j)=0.0d0 
             endif
             if (g_n(j)/g_o(j).ge.ran()) then
                ip_o(j)=ip_n(j)
                ie_o(j)=ie_n(j)
                g_o(j)=g_n(j)
                i_acc=i_acc+1
             endif
             if(g_o(j).eq.0.0d0) cycle
          
             if (iv.ge.neq.and.mod(iv,nvoid).eq.0) then
                call f_eval(nZ,p(ip_o(j)),xe(ie_o(j)),w(iw),qval,mqe,pke(ie_o(j),ip_o(j)),fg,kf,thetalept,ee,iform,f_o(j))                
                f_o(j)=f_o(j)/g_o(j)*1.e9
                r_avg=r_avg+f_o(j)
                r_err=r_err+f_o(j)**2
                i_avg=i_avg+1
             endif
          enddo
       enddo
       call addall(r_avg,r_avg_tot)
       call addall (r_err,r_err_tot)
       call addall (i_avg,i_avg_tot)
       call addall (i_acc,i_acc_tot)
       if (myrank().eq.0) then
          r_avg_tot=r_avg_tot/dble(i_avg_tot)
          r_err_tot=r_err_tot/dble(i_avg_tot)
          r_err_tot=sqrt((r_err_tot-r_avg_tot**2)/dble(i_avg_tot-1))
       !   write(6,*) 'acceptance',dble(i_acc_tot)/dble(nev*nwlk*nproc())
          write(6,*) w(iw),r_avg_tot,r_err_tot
          write(14,*) w(iw),r_avg_tot,r_err_tot
          flush(14)
       endif
    enddo
    close(14)

 !   tf=MPI_Wtime()
 !   if (myrank().eq.0) then
 !      write(6,*)'Elapsed time is',tf-ti
 !   endif
    call done() 

  end program response_ia

  subroutine f_eval(nZ,p,e,w,qval,mqe,pke,fg,kf,thetalept,ee,iform,f_o)
    use mathtool
    implicit none
    real*8, parameter :: hbarc=197.327053d0,pi= 4.0d0*atan(1.0d0)
    integer*4 :: fg,iform,nZ
    real*8 :: w,wt,qval,mqe,thetalept,ee,f_o,pke,p,kf
    real*8 :: ep,coste,pf,phi,epf,cost_te,sig,e

    f_o=0.0d0
    
    ep=sqrt(p**2+mqe**2)
    if(fg.eq.1)then
       wt=w-0.0d0
    else
       wt=w-e+mqe-ep
    endif
    cost_te=((wt+ep)**2-p**2-qval**2-mqe**2)/(2.0d0*p*qval)
    if(abs(cost_te).gt.1.0d0) then
       f_o=0.0d0
       return
    endif
    pf=sqrt(p**2+qval**2+2.0d0*qval*p*cost_te)
    epf=sqrt(mqe**2+pf**2)
    if(pf.ge.kf) then
       phi=2.0d0*pi*ran()
       call cc1(qval/hbarc,w,wt,p/hbarc,pf/hbarc,phi,ee,thetalept,iform,sig)

       f_o=p**2*pke*(dble(nZ)*sig)*epf/(p*qval)*2.0d0*pi
    endif

    return
    end subroutine


    subroutine g_eval(pj,PkE,g)
    implicit none
    real*8, parameter :: pi=acos(-1.0d0)
    real*8 :: pj,PkE,g
    g=(4.0d0*pi)*pj**2*PkE
   
    return
    end subroutine
