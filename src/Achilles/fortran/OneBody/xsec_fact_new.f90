    subroutine fortran_interface(process_id,xmn,p_4,pp_4,q_4,ff,ff_len, nspin_in,nspin_f)
      use dirac_matrices
      use libprocess_info
      use iso_c_binding
      !include "../process_info_cdef.f90"
      
      
      !implicit none
      integer :: nspin_in,nspin_f,ff_len
      real*8 :: xmn, p_4(4),pp_4(4),q_4(4)
      !real*8 :: ff1(2),ff2(2),ffa(2)
      real*8 :: ff(ff_len),ff1(2),ff2(2),ffa(2)
      complex*16,allocatable :: J_mu(:,:,:,:)

      integer*8, pointer, dimension(:) :: ids_v
      integer(c_size_t):: len,len_had
        real(c_double), pointer, dimension(:) :: masses_v
      


      type(process_info) :: process_id

      call process_id%ids(ids_v,len) ! check size of ids_v
      call process_id%masses(masses_v, len_had) ! check size of ids_v



      
      call dirac_matrices_in(xmn) !.......this has to be moved in a place where it only gets called once 
      call current_init_had(p_4,pp_4,q_4,nspin_in,nspin_f)
      call define_spinors()
      allocate(J_mu(nspin_f,nspin_in,4,2)) !...last two indices are for proton and neutron
      !....
      ff1(1)=ff(1)
      ff2(1)=ff(2)
      ffa(1)=ff(3)
      ffa(2)=0.0d0
      !....
      call det_Ja(ff1,ff2,ffa)
      call hadr_curr_matrix_el(J_mu)
      return

    end subroutine fortran_interface


! check QE_subroutine in Fortran_Interface

!..I take process info as an input: type process_info. from there I need to , library: libprocess_info. 

!... in nuclear model copy in your code: model/test--> needs to do in a module

! add call factory...and include it there. 



!...Nuclear_model_interface

!..Particle_info: contains info about particles; need to hook libpartinfo    
 


!Cmake.txt add the files there, check pke_xsec

  

