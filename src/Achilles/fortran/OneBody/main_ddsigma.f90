program response_ia
  use mathtool
  use mympi
  use dirac_matrices
  use quasi_el
  implicit none
    integer*4, parameter :: nev=200000,neq=100,nvoid=10
    integer*8, allocatable :: irn(:),irn0(:)
    real*8, parameter :: xme=0.511d0
    real*8, parameter :: small=1e-12


    integer*4 :: np0,nh0,nA,nZ,ne,nc,nw,i,j,k,nskip,iform,ie,iemax,nbox,fg,in
    integer*4, allocatable :: ip_o(:),ip_n(:),ie_o(:),ie_n(:)
    integer*4 :: nwlk
    integer*4 :: i_acc,i_avg,i_acc_tot,i_avg_tot,iv,iw
    real*8, allocatable :: w(:),wt(:),sigma(:),Q2(:)

    real*8, allocatable :: sig(:,:),sig_err(:,:)
    real*8 :: ee,thetalept,coste,eef,sig_p,sig_n,norm,mqef,dummy
    real*8 :: pmax,hp,wmax,wmin,hw,rdummy,mspec,mnuc,qval,qfm,eps,arg,mstar
    real*8 :: cost,cost_rel,dpf,kf,TA,dwfold,ffold,hwfold,hc,delta_w,espec,epf,ep
    real*8 :: r_avg,r_err,r_avg_tot,r_err_tot,mpi,xbj,he_p,he_n
    real*8, allocatable :: g_o(:),g_n(:),f_o(:),f_n(:)
    real*8 :: ti,tf
    real*8 :: phi,pf,cost_te,sin_te
    real*8 :: p_4(4),pf_4(4)
    character * 40 :: nk_fname
    character*50 :: fname,fname_pkep,fname_pken,en_char,theta_char,int_char,dec_char
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

    fname='results/C12_'//trim(en_char)//'_'//trim(theta_char)//'.out'
    fname=trim(fname)

! write cross section on file
    open(unit=14,file=fname)



    !....initialization of some stuff needs to be done before generating phase space etc    
    fname_pkep='pke12_tot.data'
    fname_pken='pke12_tot.data'
    iform=2
    call init_pke(fname_pkep,fname_pken,fg,nZ,nA,kF,iform)
    mqef=mqe/hbarc
    qfm=qval/hbarc
    thetalept=thetalept*pi/180.0d0
    coste= cos(thetalept)
    call dirac_matrices_in(mqef)
    !....end of initialization

 
    
! construct omega grid and Q2 grid
    allocate(w(nw),wt(nw),Q2(nw),sig(2,nw),sig_err(2,nw))
    wmin=0.0d0
    hw=(wmax-wmin)/dble(nw)
    do i=1,nw
       w(i)=wmin+(dble(i)-0.5d0)*hw
       eef = ee - w(i)
       Q2(i) = 2.0d0*ee*eef*(1.0d0 - coste)
    enddo

! compute the cross section in the impulse approximation
    do iw=1,nw
       qval=sqrt( Q2(iw)+ w(iw)**2 )

       do in=1,2
!
          r_avg=0.0d0
          r_err=0.0d0
          i_acc=0
          i_avg=0
          g_o=0.0d0
!..prepare importance sampling          
          do i=1,nwlk
             call setrn(irn(i))
             do while(g_o(i).le.0.0d0)
                ip_o(i)=1+int(np*ran())
                if(in.eq.1) then
                        ne=nep
                        ie_o(i)=1+int(ne*ran())
                        call g_eval(p(ip_o(i)),PkE_p(ie_o(i),ip_o(i)),g_o(i))
                elseif(in.eq.2) then
                        ne=nen
                        ie_o(i)=1+int(ne*ran())
                        call g_eval(p(ip_o(i)),PkE_n(ie_o(i),ip_o(i)),g_o(i))
                endif
             enddo
             call getrn(irn(i))
          enddo
!..start the real calculation using previously generated configurations         

          do iv=1,nev
             do j=1,nwlk
             !call setrn(irn(j))
                ip_n(j)=nint(ip_o(j)+0.05d0*np*(-1.0d0+2.0d0*ran()))
                ie_n(j)=nint(ie_o(j)+0.05d0*ne*(-1.0d0+2.0d0*ran()))
                if (ip_n(j).le.np.and.ip_n(j).ge.1.and.ie_n(j).le.ne.and.ie_n(j).ge.1) then
                   if(in.eq.1) call g_eval(p(ip_n(j)),PkE_p(ie_n(j),ip_n(j)),g_n(j))
                   if(in.eq.2) call g_eval(p(ip_n(j)),PkE_n(ie_n(j),ip_n(j)),g_n(j))
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

                   !...we generate the rest of the phasespace
                   phi=2.0d0*pi*ran()
                   cost_te=-1.0d0+2.0d0*ran()
                   sin_te=sqrt(1.0d0-cost_te**2)
                   ep=sqrt(p(ip_o(j))**2+mqe**2)
                   p_4(1)=ep
                   p_4(2)=p(ip_o(j))*sin_te*cos(phi)
                   p_4(3)=p(ip_o(j))*sin_te*sin(phi)
                   p_4(4)=p(ip_o(j))*cost_te
                   pf=sqrt(p(ip_o(j))**2+qval**2+2.0d0*qval*p(ip_o(j))*cost_te)
                   !if(pf.lt.kf) cycle
                   epf=sqrt(mqe**2+pf**2)
                   pf_4(1)=epf
                   pf_4(2)=p(ip_o(j))*sin_te*cos(phi)
                   pf_4(3)=p(ip_o(j))*sin_te*sin(phi)
                   pf_4(4)=p(ip_o(j))*cost_te+qval
                  
                   call f_eval(in,p_4,pf_4,ie_o(j),ip_o(j),w(iw),qval,thetalept,ee,f_o(j))

                   f_o(j)=f_o(j)/g_o(j)*1.e6
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
             sig(in,iw)=r_avg_tot
             sig_err(in,iw)=r_err_tot
          endif
       enddo
       if (myrank().eq.0) then
         xbj=Q2(iw)/(2.0d0*mp*w(iw))
         !   write(6,*) 'acceptance',dble(i_acc_tot)/dble(nev*nwlk*nproc())
         write(6,*) w(iw),xbj,sum(sig(:,iw)),sum(sig_err(:,iw))
         write(14,*) w(iw),xbj,sum(sig(:,iw)),sig(1,iw),sig(2,iw)
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



    subroutine g_eval(pj,PkE,g)
    implicit none
    real*8, parameter :: pi=acos(-1.0d0)
    real*8 :: pj,PkE,g
    g=(4.0d0*pi)*pj**2*PkE
   
    return
    end subroutine

