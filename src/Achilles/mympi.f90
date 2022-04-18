!$Id: mympi.f90,v 1.6 2013/12/10 21:21:53 nuclear Exp $
module mympi
!
! load balancing mpi routines -- these are rewrites by K.E. Schmidt
! of the mpi routines written by Michael A. Lee and I. Lomonosov for the
! parallel version of the Schmidt and Lee electronic structure GFMC code
! 
   implicit none
   include 'mpif.h'
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: irank,iproc,npart
   integer(kind=i4), private, save :: ncheck,nbal,totalnw
   integer(kind=i4), private, save, allocatable :: numfig(:)
   integer(kind=i4), private, save, allocatable :: nsend(:),nrecv(:),nmoves(:)
   integer(kind=i4), private, save, allocatable :: nwmax(:),nwmin(:)
   real(kind=r8),    private, save :: dobal
   real(kind=r8),    private, save :: t0,t1,tbal
   logical,          private, save :: ivmc

   interface bcast  ! broadcast from process 0
      module procedure bcasti1,bcasti81,bcasti1d,bcasti2d,bcasti3d,bcasti4d
      module procedure bcastr1,bcastr1d,bcastr2d,bcastr3d,bcastr4d
      module procedure bcastc1,bcastc1d,bcastc2d,bcastc3d,bcastc4d
      module procedure bcastb,bcastb1d
      module procedure bcastchar
   end interface bcast

   interface addall  ! return sum to process 0
      module procedure addalli1,addalli1d
      module procedure addallr1,addallr1d,addallr2d
      module procedure addallc1,addallc1d,addallc2d
   end interface addall

   interface gather  ! gather to process 0
      module procedure gatheri1,gatheri1d
      module procedure gatherr1,gatherr1d
      module procedure gatherc1,gatherc1d
   end interface gather

   interface minall  ! allreduce with mpi_min
      module procedure minallr1
   end interface minall

contains

subroutine init0()  ! call this before anything else
   integer(kind=i4) :: ierror
   call mpi_init(ierror)
   call mpi_comm_rank(mpi_comm_world,irank,ierror)
   call mpi_comm_size(mpi_comm_world,iproc,ierror)
   allocate(numfig(0:iproc-1))
   if (mpi_integer8.eq.0) call stop('mpi_integer8 not defined ',0)
   allocate(nsend(0:iproc-1),nrecv(0:iproc-1),nmoves(0:iproc-1))
   allocate(nwmax(0:iproc-1),nwmin(0:iproc-1))
   nsend=0
   nrecv=0
   nmoves=0
   nwmax=0
   nwmin=0
   ncheck=0
   nbal=0
   totalnw=0
   call cpu_time(t0)
end subroutine init0

subroutine init1(npartin,dobalin,ivmcin)  ! call this when everyone knows these
   integer(kind=i4) :: npartin
   real(kind=r8) :: dobalin
   logical :: ivmcin
   npart=npartin
   dobal=dobalin
   ivmc=ivmcin
end subroutine init1

subroutine done()  ! wrapper for finalize routine
   integer(kind=i4) :: ierror
   call mpi_finalize(ierror)
end subroutine done

subroutine bcasti1(i)
   integer(kind=i4) :: i,ierror
   call mpi_bcast(i,1,mpi_integer,0,mpi_comm_world,ierror)
   return
end subroutine bcasti1

subroutine bcasti81(i)
   integer(kind=i8) :: i
   integer(kind=i4) :: ierror
   call mpi_bcast(i,1,mpi_integer8,0,mpi_comm_world,ierror)
   return
end subroutine bcasti81

subroutine bcasti1d(i)
   integer(kind=i4) :: i(:),ierror
   call mpi_bcast(i,size(i),mpi_integer,0,mpi_comm_world,ierror)
   return
end subroutine bcasti1d

subroutine bcasti2d(i)
   integer(kind=i4) :: i(:,:),ierror
   call mpi_bcast(i,size(i),mpi_integer,0,mpi_comm_world,ierror)
   return
end subroutine bcasti2d

subroutine bcasti3d(i)
   integer(kind=i4) :: i(:,:,:),ierror
   call mpi_bcast(i,size(i),mpi_integer,0,mpi_comm_world,ierror)
   return
end subroutine bcasti3d

subroutine bcasti4d(i)
   integer(kind=i4) :: i(:,:,:,:),ierror
   call mpi_bcast(i,size(i),mpi_integer,0,mpi_comm_world,ierror)
   return
end subroutine bcasti4d

subroutine bcastr1(r)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r
   call mpi_bcast(r,1,mpi_double_precision,0,mpi_comm_world,ierror)
   return
end subroutine bcastr1

subroutine bcastr1d(r)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r(:)
   call mpi_bcast(r,size(r),mpi_double_precision,0,mpi_comm_world,ierror)
   return
end subroutine bcastr1d

subroutine bcastr2d(r)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r(:,:)
   call mpi_bcast(r,size(r),mpi_double_precision,0,mpi_comm_world,ierror)
   return
end subroutine bcastr2d

subroutine bcastr3d(r)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r(:,:,:)
   call mpi_bcast(r,size(r),mpi_double_precision,0,mpi_comm_world,ierror)
   return
end subroutine bcastr3d

subroutine bcastr4d(r)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r(:,:,:,:)
   call mpi_bcast(r,size(r),mpi_double_precision,0,mpi_comm_world,ierror)
   return
end subroutine bcastr4d

subroutine bcastc1(c)
   integer(kind=i4) :: ierror
   complex(kind=r8) :: c
   call mpi_bcast(c,1,mpi_double_precision,0,mpi_comm_world,ierror)
   return
end subroutine bcastc1

subroutine bcastc1d(c)
   integer(kind=i4) :: ierror
   complex(kind=r8) :: c(:)
   call mpi_bcast(c,size(c),mpi_double_complex,0,mpi_comm_world,ierror)
   return
end subroutine bcastc1d

subroutine bcastc2d(c)
   integer(kind=i4) :: ierror
   complex(kind=r8) :: c(:,:)
   call mpi_bcast(c,size(c),mpi_double_complex,0,mpi_comm_world,ierror)
   return
end subroutine bcastc2d

subroutine bcastc3d(c)
   integer(kind=i4) :: ierror
   complex(kind=r8) :: c(:,:,:)
   call mpi_bcast(c,size(c),mpi_double_complex,0,mpi_comm_world,ierror)
   return
end subroutine bcastc3d

subroutine bcastc4d(c)
   integer(kind=i4) :: ierror
   complex(kind=r8) :: c(:,:,:,:)
   call mpi_bcast(c,size(c),mpi_double_complex,0,mpi_comm_world,ierror)
   return
end subroutine bcastc4d

subroutine bcastb(b)
   integer(kind=i4) :: ierror
   logical :: b
   call mpi_bcast(b,1,mpi_logical,0,mpi_comm_world,ierror)
   return
end subroutine bcastb

subroutine bcastb1d(b)
   integer(kind=i4) :: ierror
   logical :: b(:)
   call mpi_bcast(b,size(b),mpi_logical,0,mpi_comm_world,ierror)
   return
end subroutine bcastb1d

subroutine bcastchar(c)
   integer(kind=i4) :: ierror
   character(len=*) :: c
   call mpi_bcast(c,len(c),mpi_character,0,mpi_comm_world,ierror)
   return
end subroutine bcastchar

function myrank()  ! which process am I?
   integer(kind=i4) :: myrank
   myrank=irank
end function myrank

function nproc()  ! How many of use are there anyway?
   integer(kind=i4) :: nproc
   nproc=iproc
end function nproc

subroutine barrier()  ! wrapper for mpi_barrier
   integer(kind=i4) :: ierror
   call mpi_barrier(mpi_comm_world,ierror)
end subroutine barrier

subroutine addalli1(i,isum)
   integer(kind=i4) :: ierror,i,isum
   call mpi_reduce(i,isum,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierror)
   return
end subroutine addalli1

subroutine addalli1d(i,isum)
   integer(kind=i4) :: ierror,i(:),isum(:)
   call mpi_reduce(i,isum,size(i),mpi_integer,mpi_sum,0,mpi_comm_world,ierror)
   return
end subroutine addalli1d

subroutine addallr1(r,rsum)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r,rsum
   call mpi_reduce(r,rsum,1,mpi_double_precision,mpi_sum,0, &
      mpi_comm_world,ierror)
   return
end subroutine addallr1

subroutine addallr1d(r,rsum)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r(:),rsum(:)
   call mpi_reduce(r,rsum,size(r),mpi_double_precision,mpi_sum,0, &
      mpi_comm_world,ierror)
   return
end subroutine addallr1d

subroutine addallr2d(r,rsum)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r(:,:),rsum(:,:)
   call mpi_reduce(r,rsum,size(r),mpi_double_precision,mpi_sum,0, &
      mpi_comm_world,ierror)
   return
end subroutine addallr2d

subroutine addallc1(c,csum)
   integer(kind=i4) :: ierror
   complex(kind=r8) :: c,csum
   call mpi_reduce(c,csum,1,mpi_double_complex,mpi_sum,0, &
      mpi_comm_world,ierror)
   return
end subroutine addallc1

subroutine addallc1d(c,csum)
   integer(kind=i4) :: ierror
   complex(kind=r8) :: c(:),csum(:)
   call mpi_reduce(c,csum,size(c),mpi_double_complex,mpi_sum,0, &
      mpi_comm_world,ierror)
   return
end subroutine addallc1d

subroutine addallc2d(c,csum)
   integer(kind=i4) :: ierror
   complex(kind=r8) :: c(:,:),csum(:,:)
   call mpi_reduce(c,csum,size(c),mpi_double_complex,mpi_sum,0, &
      mpi_comm_world,ierror)
   return
end subroutine addallc2d

subroutine gatheri1(i,igather)
   integer(kind=i4) :: i,igather(:),ierror
   call mpi_gather(i,1,mpi_integer,igather,1,mpi_integer,0, &
      mpi_comm_world,ierror)
   return
end subroutine gatheri1

subroutine gatheri1d(i,igather)
   integer(kind=i4) :: i(:),igather(:,:),ierror
   call mpi_gather(i,size(i),mpi_integer,igather,size(i),mpi_integer,0, &
      mpi_comm_world,ierror)
   return
end subroutine gatheri1d

subroutine gatherr1(r,rgather)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r,rgather(:)
   call mpi_gather(r,1,mpi_double_precision,rgather,1, &
      mpi_double_precision,0,mpi_comm_world,ierror)
   return
end subroutine gatherr1

subroutine gatherr1d(r,rgather)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r(:),rgather(:,:)
   call mpi_gather(r,size(r),mpi_double_precision,rgather,size(r) &
      ,mpi_double_precision,0,mpi_comm_world,ierror)
   return
end subroutine gatherr1d

subroutine gatherc1(c,cgather)
   integer(kind=i4) :: ierror
   complex(kind=r8) :: c,cgather(:)
   call mpi_gather(c,1,mpi_double_complex,cgather,1, &
      mpi_double_complex,0,mpi_comm_world,ierror)
   return
end subroutine gatherc1

subroutine gatherc1d(c,cgather)
   integer(kind=i4) :: ierror
   complex(kind=r8) :: c(:),cgather(:,:)
   call mpi_gather(c,size(c),mpi_double_complex,cgather,size(c) &
      ,mpi_double_complex,0,mpi_comm_world,ierror)
   return
end subroutine gatherc1d

subroutine minallr1(r,rmin)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r,rmin
   call mpi_allreduce(r,rmin,1,mpi_double_precision,mpi_min,mpi_comm_world,ierror)
   return
end subroutine minallr1
end module mympi

subroutine stop(msg,icode)
   use mympi
   implicit none
   integer, parameter :: i4=selected_int_kind(9)
   integer(kind=i4) :: icode,ierror
   character*(*) :: msg
   write ( 6, '(/,1x,78(''*''),/,'' * '',''Stopping because of an error'', &
      t79,''*'',/,'' * '',a,i10,t79,''*'',/,1x,78(''*''),/)' ) msg,icode
   call mpi_abort(mpi_comm_world,ierror)
   stop 9999
   return
end subroutine stop

