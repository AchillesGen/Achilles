module mathtool
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i8), private, save :: irn
   
contains

subroutine locale(v,nv,val,iv)
   implicit none
   integer*4 :: i,nv,iv
   real*8 :: v(nv),val
   iv=1
   do i=1,nv
    if(abs(v(i)-val).lt.abs(v(i)-v(iv))) iv=i
   enddo
   return
end subroutine

subroutine interpolint(xin,yin,nin,x,y,np)
   integer*4 :: nin,np,npp,i,j
   real*8 :: xin(nin),yin(nin),x,y,dy
  
   call locate(xin,nin,x,j)
   if (j.eq.0) j=1
   npp=min(np,iabs(nin-j))
   if (npp.eq.nin) y=yin(nin)
   call polint(xin(j:j+npp),yin(j:j+npp),npp+1,x,y,dy)
end subroutine interpolint

!===================================================================================
!Given arrays xa and ya, each of length n, and given a value x, this routine returns a value y,
! and an error estimate dy. If P (x) is the polynomial of degree N − 1 such that 
!P(xai) = yai,i = 1,...,n, then the returned value y = P(x).
!===================================================================================
subroutine polint(xa,ya,n,x,y,dy)
   integer*4,parameter :: NMAX=10
   integer*4 :: n,i,m,ns
   real*8 :: dy,x,y,xa(n),ya(n)
   real*8 :: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
   ns=1
   dif=abs(x-xa(1))
   do i=1,n
      dift=abs(x-xa(i))
      if (dift.lt.dif) then
         ns=i
         dif=dift
      endif
      c(i)=ya(i)
      d(i)=ya(i)
   enddo
   y=ya(ns)
   ns=ns-1
   do m=1,n-1
      do i=1,n-m
         ho=xa(i)-x
         hp=xa(i+m)-x
         w=c(i+1)-d(i)
         den=ho-hp
         if(den.eq.0.)stop 'failure in polint'
         den=w/den
         d(i)=hp*den
         c(i)=ho*den
      enddo
      if (2*ns.lt.n-m)then
         dy=c(ns+1)
      else
         dy=d(ns)
         ns=ns-1
      endif
      y=y+dy
   enddo
   return
end subroutine polint

!===================================================================================
!Given an array xx(1:n), and given a value x, returns a value j such that x is between xx(j)
! and xx(j+1). xx(1:n) must be monotonic, either increasing or decreasing. j=0 or j=n is 
!returned to indicate that x is out of range. 
!===================================================================================
subroutine locate(xx,n,x,j) 
   integer*4 :: j,n,jl,jm,ju
   real*8 :: x,xx(n)

   jl=0 
   ju=n+1
   do while (ju-jl.gt.1)  
      jm=(ju+jl)/2
      if ((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm))) then
         jl=jm 
      else
         ju=jm 
      endif
   enddo
   if (x.eq.xx(1))then
      j=1
   else if(x.eq.xx(n))then
      j=n-1 
   else
      j=jl 
   endif
   return
end subroutine locate


!===================================================================================
!Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f (xi), with
!x1 < x2 < . . . < xN , and given values yp1 and ypn for the first derivative of the inter-
!polating function at points 1 and n, respectively, this routine returns an array y2(1:n) of
!length n which contains the second derivatives of the interpolating function at the tabulated
!points xi . If nat1 and/or natn are equal to 1, the routine is signaled to set
!the corresponding boundary condition for a natural spline, with zero second derivative on
!that boundary. If nat1 and/or natn are equal to 2, the routine automatically computes the
!first derivative at the boundaries with the method of five points. 
!===================================================================================
subroutine spline_init_nr(x,y,n,yp1,ypn,y2,nat1,natn)

   implicit none

   integer*4, parameter :: NMAX=15000        !Parameter: NMAX is the largest anticipated value of n.
   integer*4 :: i,k,n,nat1,natn                      
   real*8 :: yp1,ypn,fp1,fpn,x(n),y(n),y2(n)
   real*8 :: p,qn,sig,un,u(NMAX)
 
!The lower boundary condition is set either to be
!“natural” (nat1=1) or else to have a specified first derivative.
   if (nat1.eq.1) then
    y2(1)=0.0d0
    u(1)=0.0d0
   else
!If required (nat1=2), the first derivative at the left border is computed with five points
    if (nat1.eq.2) then
     fp1=(-3.0d0*y(5)+16.0d0*y(4)-36.0d0*y(3)+48d0*y(2)-25.0d0*y(1))/12.0d0/(x(2)-x(1))
    else
     fp1=yp1
    endif
    y2(1)=-0.5d0
    u(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-fp1)
   endif
!This is the decomposition loop of the tridiagonal
!algorithm. y2 and u are used for temporary
!storage of the decomposed factors.
   do i=2,n-1
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2.0d0
    y2(i)=(sig-1.)/p
    u(i)=(6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))&
  & /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
   enddo

!The upper boundary condition is set either to be
!“natural” (natn=1), or else to have a specified first derivative.
   if (natn.eq.1) then
    qn=0.0d0
    un=0.0d0
   else
!If required (natn=2), the first derivative at the right border is computed with five points
    if(natn.eq.2) then
     fpn=(25.0d0*y(n)-48.0d0*y(n-1)+36.0d0*y(n-2)-16.0d0*y(n-3)+3.0d0*y(n-4))/12.0d0/(x(2)-x(1))
    else
     fpn=ypn
    endif
    qn=0.5d0
    un=(3.0d0/(x(n)-x(n-1)))*(fpn-(y(n)-y(n-1))/(x(n)-x(n-1)))
   endif
   y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)

!This is the backsubstitution loop of the tridiago-
!nal algorithm.
   do k=n-1,1,-1
    y2(k)=y2(k)*y2(k+1)+u(k)
   enddo 

   return
end subroutine spline_init_nr


!Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
!xai ’s in order), and given the array y2a(1:n), which is the output from spline above,
!and given a value of x, this routine returns a cubic-spline interpolated value y.
subroutine splint_nr(xa,ya,y2a,n,x,y)
   integer*4 :: n,k,khi,klo
   real*8 :: x,y,xa(n),y2a(n),ya(n),a,b,h

!We will find the right place in the table by means of bisection.    
   klo=1
   khi=n
!This is optimal if sequential calls to this routine are at random
!values of x. If sequential calls are in order, and closely
!spaced, one would do better to store previous values of
!klo and khi and test if they remain appropriate on the next call.
1  if (khi-klo.gt.1) then
    k=(khi+klo)/2
    if(xa(k).gt.x)then
     khi=k
    else
     klo=k
    endif
    goto 1
   endif
!klo and khi now bracket the input value of x.
   h=xa(khi)-xa(klo)
   if (h.eq.0.) then
    write(*,*) 'bad xa input in splint' ! The xa’s must be distinct.
   endif
   a=(xa(khi)-x)/h
!Cubic spline polynomial is now evaluated.
   b=(x-xa(klo))/h
   y=a*ya(klo)+b*ya(khi)+&
   &((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0d0
   return
end subroutine splint_nr

!===================================================================================
! questa subroutine calcola la derivata della funzione f e restituisce fp
!===================================================================================
subroutine der_3p(xr,nr,nstep,f,fp,fpp)
   integer(kind=i4) :: i,nr,nstep
   real(kind=r8) :: hr,xr(nr),f(nr),fp(nr),fpp(nr)
   hr=xr(2)-xr(1)    
! 1st derivative
   fp(1)=(-1.0_r8*f(3)+4.0_r8*f(2)-3.0_r8*f(1))/(2.0_r8*hr)
   do i=2,nstep-1
      fp(i)=(f(i+1)-f(i-1))/(2.0_r8*hr)
   enddo
   fp(nstep)=(3.0_r8*f(nstep)-4.0_r8*f(nstep-1)+1.0_r8*f(nstep-2))/(2.0_r8*hr)
! 2nd derivative
   fpp(1)=(f(3)-2.0_r8*f(2)+f(1))/hr**2
   do i=2,nstep-1
      fpp(i)=(f(i+1)-2.0_r8*f(i)+f(i-1))/hr**2
   enddo
   fpp(nstep)=(f(nstep)-2.0_r8*f(nstep-1)+f(nstep-2))/hr**2
   return
end subroutine der_3p

!===================================================================================
! questa subroutine calcola la derivata della funzione f e restituisce fp
! usando cinque punti per differenziare
!===================================================================================
subroutine der_5p(xr,nr,nstep,f,fp,fpp)
   integer(kind=i4) :: i,nr,nstep
   real(kind=r8) :: hr,xr(nr),f(nr),fp(nr),fpp(nr)
   hr=xr(2)-xr(1)    
! 1st derivative
   fp(1)=(-3.0_r8*f(5)+16.0_r8*f(4)-36.0_r8*f(3)+48.0_r8*f(2)-25.0_r8*f(1))/(12.0_r8*hr)
   fp(2)=( 1.0_r8*f(5)- 6.0_r8*f(4)+18.0_r8*f(3)-10.0_r8*f(2)- 3.0_r8*f(1))/(12.0_r8*hr)
   do i=3,nstep-2
      fp(i)=(-f(i+2)+8.0_r8*f(i+1)-8.0_r8*f(i-1)+f(i-2))/(12.0_r8*hr)
   enddo
   fp(nstep-1)=( 3.0_r8*f(nstep)+10.0_r8*f(nstep-1)-18.0_r8*f(nstep-2)+ 6.0_r8*f(nstep-3)-&
   1.0_r8*f(nstep-4))/(12.0_r8*hr)
   fp(nstep)  =(25.0_r8*f(nstep)-48.0_r8*f(nstep-1)+36.0_r8*f(nstep-2)-16.0_r8*f(nstep-3)+&
   3.0_r8*f(nstep-4))/(12.0_r8*hr)
! 2nd derivative
   fpp(1)=( 11.0_r8*f(5)-56.0_r8*f(4)+114.0_r8*f(3)-104.0_r8*f(2)+35.0_r8*f(1))/(12.0_r8*hr**2)
   fpp(2)=(- 1.0_r8*f(5)+ 4.0_r8*f(4)+  6.0_r8*f(3)- 20.0_r8*f(2)+11.0_r8*f(1))/(12.0_r8*hr**2)
   do i=3,nstep-2
      fpp(i)=(-f(i+2)+16.0_r8*f(i+1)-30.0_r8*f(i)+16.0_r8*f(i-1)-f(i-2))/(12.0_r8*hr**2)
   enddo
   fpp(nstep-1)=(11.0_r8*f(nstep)- 20.0_r8*f(nstep-1)+  6.0_r8*f(nstep-2)+ 4.0_r8*f(nstep-3)-&
   1.0_r8*f(nstep-4))/(12.0_r8*hr**2)
   fpp(nstep)  =(35.0_r8*f(nstep)-104.0_r8*f(nstep-1)+114.0_r8*f(nstep-2)-56.0_r8*f(nstep-3)+&
   11.0_r8*f(nstep-4))/(12.0_r8*hr**2)
   return
end subroutine der_5p

subroutine der_7p(xr,nr,nstep,f,fp,fpp)
   integer(kind=i4) :: i,nr,nstep
   real(kind=r8) :: hr,xr(nr),f(nr),fp(nr),fpp(nr)
   hr=xr(2)-xr(1)
! 1st derivative
   fp(1)=(-10.0_r8*f(7)+72.0_r8*f(6)-225.0_r8*f(5)+400.0_r8*f(4)-450.0_r8*f(3)+360.0_r8*f(2)-&
   147.0_r8*f(1))/(60.0_r8*hr)
   fp(2)=(  2.0_r8*f(7)-15.0_r8*f(6)+ 50.0_r8*f(5)-100.0_r8*f(4)+150.0_r8*f(3)- 77.0_r8*f(2)-&
   10.0_r8*f(1))/(60.0_r8*hr)
   fp(3)=(- 1.0_r8*f(7)+ 8.0_r8*f(6)- 30.0_r8*f(5)+ 80.0_r8*f(4)- 35.0_r8*f(3)- 24.0_r8*f(2)+&
   2.0_r8*f(1))/(60.0_r8*hr)
   do i=4,nstep-3
      fp(i)=(f(i+3)-9.0_r8*f(i+2)+45.0_r8*f(i+1)-45.0_r8*f(i-1)+9.0_r8*f(i-2)-f(i-3))/(60.0_r8*hr)
   enddo
   fp(nstep-2)=(-  2.0_r8*f(nstep)+ 24.0_r8*f(nstep-1)+ 35.0_r8*f(nstep-2)- 80.0_r8*f(nstep-3)+&
   30.0_r8*f(nstep-4)-  8.0_r8*f(nstep-5)+  1.0_r8*f(nstep-6))/(60.0_r8*hr)
   fp(nstep-1)=(  10.0_r8*f(nstep)+ 77.0_r8*f(nstep-1)-150.0_r8*f(nstep-2)+100.0_r8*f(nstep-3)-&
   50.0_r8*f(nstep-4)+ 15.0_r8*f(nstep-5)-  2.0_r8*f(nstep-6))/(60.0_r8*hr)
   fp(nstep)  =( 147.0_r8*f(nstep)-360.0_r8*f(nstep-1)+450.0_r8*f(nstep-2)-400.0_r8*f(nstep-3)+&
   225.0_r8*f(nstep-4)- 72.0_r8*f(nstep-5)+ 10.0_r8*f(nstep-6))/(60.0_r8*hr)
! 2nd derivative
   fpp(1)=( 137.0_r8*f(7)-972.0_r8*f(6)+2970.0_r8*f(5)-5080.0_r8*f(4)+5265.0_r8*f(3)-&
   3132.0_r8*f(2)+812.0_r8*f(1))/(180.0_r8*hr**2)
   fpp(2)=(- 13.0_r8*f(7)+ 93.0_r8*f(6)- 285.0_r8*f(5)+ 470.0_r8*f(4)- 255.0_r8*f(3)-&
   147.0_r8*f(2)+137.0_r8*f(1))/(180.0_r8*hr**2)
   fpp(3)=(   2.0_r8*f(7)- 12.0_r8*f(6)+  15.0_r8*f(5)+ 200.0_r8*f(4)- 420.0_r8*f(3)+&
   228.0_r8*f(2)- 13.0_r8*f(1))/(180.0_r8*hr**2)
   do i=4,nr-3
      fpp(i)=(2.0_r8*f(i+3)-27.0_r8*f(i+2)+270.0_r8*f(i+1)-490_r8*f(i)+270.0_r8*f(i-1)-&
   27.0_r8*f(i-2)+2.0_r8*f(i-3))/(180.0_r8*hr**2)
   enddo
   fpp(nstep-2)=(- 13.0_r8*f(nstep)+ 228.0_r8*f(nstep-1)- 420.0_r8*f(nstep-2)+&
   200.0_r8*f(nstep-3)+  15.0_r8*f(nstep-4)- 12.0_r8*f(nstep-5)+  2.0_r8*f(nstep-6))/(180.0_r8*hr**2)
   fpp(nstep-1)=( 137.0_r8*f(nstep)- 147.0_r8*f(nstep-1)- 255.0_r8*f(nstep-2)+&
   470.0_r8*f(nstep-3)- 285.0_r8*f(nstep-4)+ 93.0_r8*f(nstep-5)- 13.0_r8*f(nstep-6))/(180.0_r8*hr**2)
   fpp(nstep)  =( 812.0_r8*f(nstep)-3132.0_r8*f(nstep-1)+5265.0_r8*f(nstep-2)-&
   5080.0_r8*f(nstep-3)+2970.0_r8*f(nstep-4)-972.0_r8*f(nstep-5)+137.0_r8*f(nstep-6))/(180.0_r8*hr**2)
   return
end subroutine der_7p


!===================================================================================
! questa subroutine fornisce un'approssimazione della delta di dirac
! con una gaussiana
!===================================================================================
function fdelta( x, sigma )
   implicit none
   real*8 :: pi, cnorm, fdelta, sigma, x

   pi=acos(-1.0d0)
   cnorm=1.0d0/sqrt(pi)/sigma
   fdelta=cnorm*exp(-(x/sigma)**2)

   return
end function fdelta


!===================================================================================
! Given a function func defined on the interval from x1-x2 subdivide the interval 
! into n equally spaced segments, and search for zero crossings of the function. 
! nb is returned as the number of bracketing pairs xb1(1:nb), xb2(1:nb) that are 
! found. xb1 and xb2 are pointers to arrays of length nb that are dynamically 
! allocated by the routine.
!===================================================================================
subroutine brak_ale(func,fd,nd,x1,x2,n,xb1,xb2,nb)
   implicit none
   integer*4 :: nd,n,nb
   real*8 :: fd(nd),x1,x2,xb1(n),xb2(n)
   integer*4 :: i
   real*8 :: x,dx,fc,fp
   interface 
   function func(x,fd,nd)
       implicit none
       integer, intent (in) :: nd
       real*8, intent(in) :: x,fd(nd)  
       real*8 :: func
   end function func 
   end interface

   nb=0
   dx=(x2-x1)/dble(n)
   x=x1
   fp=func(x,fd,nd)
   do i=1,n
      x=x+dx
      fc=func(x,fd,nd)
      if (fc*fp.le.0.0d0) then 
         nb=nb+1
         xb1(nb)=x-dx
         xb2(nb)=x
      endif
      fp=fc
   enddo

   return      
end subroutine

!===================================================================================
! Using bisection, find the root of a function func known to lie between x1 and x2. 
! The root, returned as rtbis, will be refined until its accuracy is ±xacc.
! Parameter: MAXIT is the maximum allowed number of bisections.
!===================================================================================
subroutine bisec_ale(rtbis,func,fd,nd,x1,x2,xacc)
   integer*4 :: nd,jmax
   real*8 :: rtbis,fd(nd),x1,x2,xacc
   interface 
   function func(x,fd,nd)
       implicit none
       integer, intent (in) :: nd
       real*8, intent(in) :: x,fd(nd)  
       real*8 :: func
   end function func 
   end interface

   integer*4 :: j
   real*8 :: dx,f,fmid,xmid

   jmax=100 ! Maximum allowed number of bisections

   fmid=func(x2,fd,nd)
   f=func(x1,fd,nd)
   if (f*fmid.ge.0.0d0) then
      write(6,*) 'Error in bisec_ale: root must be bracketed'
      stop
   endif
   if (f.lt.0.0d0) then ! Orient the search so that f>0 lies at x+dx.
      rtbis=x1
      dx=x2-x1
   else
      rtbis=x2
      dx=x1-x2
   endif
   do j=1,jmax ! Bisection loop.
      dx=dx*0.5d0
      xmid=rtbis+dx
      fmid=func(xmid,fd,nd)
      if (fmid.le.0.0d0) rtbis=xmid
      if (abs(dx).lt.xacc.or.fmid.eq.0.0d0) return
   enddo
   write(6,*) 'Error in bisec_ale: not enough bisections'
   
   return

end subroutine bisec_ale

subroutine setrn(irnin)
    integer, parameter :: i8=selected_int_kind(15)
    integer, parameter :: r8=selected_real_kind(15,9)
    integer(kind=i8) :: irnin
    irn=irnin
end subroutine setrn

function rgauss()
    implicit none
    integer, parameter :: i8=selected_int_kind(15)
    integer, parameter :: r8=selected_real_kind(15,9)
    real(kind=r8) :: pi
    real(kind=r8) :: rgauss,x1,x2
    x1=ran()
    x2=ran()
    pi=4.0_r8*atan(1.0_r8)
    rgauss=sqrt(-2.0_r8*log(x1))*cos(2.0_r8*pi*x2)
    return
end function rgauss
    
function ran()
    implicit none
    integer, parameter :: i4=selected_int_kind(9)
    integer, parameter :: i8=selected_int_kind(15)
    integer, parameter :: r8=selected_real_kind(15,9)
    integer(kind=i8),  parameter :: mask24 = ishft(1_i8,24)-1
    integer(kind=i8),  parameter :: mask48 = ishft(1_i8,48_i8)-1_i8
    real(kind=r8),  parameter :: twom48=2.0_r8**(-48)
    integer(kind=i8),  parameter :: mult1 = 44485709377909_i8
    integer(kind=i8),  parameter :: m11 = iand(mult1,mask24)
    integer(kind=i8),  parameter :: m12 = iand(ishft(mult1,-24),mask24)
    integer(kind=i8),  parameter :: iadd1 = 96309754297_i8
    real(kind=r8) :: ran
    integer(kind=i8) :: is1,is2
    is2 = iand(ishft(irn,-24),mask24)
    is1 = iand(irn,mask24)
    irn = iand(ishft(iand(is1*m12+is2*m11,mask24),24)+is1*m11+iadd1,mask48)
    ran = ior(irn,1_i8)*twom48
    return
end function ran

subroutine getrn(irnout)
    integer, parameter :: i8=selected_int_kind(15)
    integer, parameter :: r8=selected_real_kind(15,9)
    integer(kind=i8) :: irnout
    irnout=irn
end subroutine getrn 

end module mathtool
    
