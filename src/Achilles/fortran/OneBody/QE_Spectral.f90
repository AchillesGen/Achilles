module QE_Spectral_Function
    use nuclear_model   
    use libutilities    
    implicit none

   type, extends(model) :: qe_sf_model
       contains



    function init_model(name) bind(C, name="pke_files.txt")
        use iso_c_binding
        use libutilities
        implicit none

        type(c_ptr), intent(in), value :: name
        logical :: init_model
        character(len=:), allocatable :: fname

        fname = c2fstring(name)

        init_model = model_ptr%init(fname)

        call init(constants) ! load constants
        call dirac_matrices_in(constants%mqe) 


    end function



!    function test_init(self, filename)

!        class(test), intent(inout) :: self
!        character(len=*), intent(in) :: filename
!        logical :: test_init
!        test_init = .true.
!        call init(constants) ! load constants
!        call dirac_matrices_in(constants%mqe)      ! not sure where I should include this   
!    end function


! NuclearModel.cc




    subroutine currents(pids_in, pids_out, moms, nin, nout, qvec, ff, len_ff, cur, nspin, nlorentz) bind(C, name="GetCurrents")
        use iso_c_binding
        use libvectors
        use dirac_matrices

        implicit none

        real(c_double), intent(in), dimension(4), target :: qvec
        integer(c_size_t), intent(in), value :: len_ff, nin, nout, nspin, nlorentz
        complex(c_double_complex), dimension(nlorentz, nspin), intent(out) :: cur
        !complex(c_double_complex), dimension(len_ff), intent(in), target :: ff !...I am changing this, they are real
        real(c_double), dimension(len_ff), intent(in), target :: ff
        


        ! C++ is row major, while Fortran is column major
        ! This means the moms passed in are transposed
        integer(c_int), dimension(nin), intent(in) :: pids_in
        integer(c_int), dimension(nout), intent(in) :: pids_out
        real(c_double), intent(in), dimension(4, nin+nout) :: moms
        type(fourvector), dimension(nin) :: mom_in
        type(fourvector), dimension(nout) :: mom_out
        type(fourvector) :: qvector
        integer(c_size_t) :: i
        !
        real(c_double), dimension(2), intent(in), target :: ff1,ff2,ffa
        complex(c_double_complex), dimension(2,2, nlorentz) :: J_mu
        

        do i=1,nin
            mom_in(i) = fourvector(moms(1, i), moms(2, i), moms(3, i), moms(4, i))
        enddo

        do i=1,nout
            mom_out(i) = fourvector(moms(1, nin+i), moms(2, nin+i), moms(3, nin+i), moms(4, nin+i))
        enddo

        qvector = fourvector(qvec(1), qvec(2), qvec(3), qvec(4))
        call current_init_had(mom_in,mom_out,qvec)
        
        call define_spinors()
        ff1(1)=ff(1)
        ff2(1)=ff(2)
        ffa(1)=ff(3)
        ffa(2)=0.0d0

        call det_Ja(ff1,ff2,ffa)
        call hadr_curr_matrix_el(J_mu)

        do i=1,2
           do j=1,2
              cur(i*j,:)= J_mu(j,i,:)
            enddo   
        enddo
     return
    end subroutine






    subroutine test_currents(self, pids_in, mom_in, nin, pids_out, mom_out, nout, qvec, ff, len_ff, cur, nspin, nlorentz)
        use iso_c_binding
        use libvectors
        use dirac_matrices
        class(test), intent(inout) :: self
        integer(c_size_t), intent(in), value :: nin, nout, len_ff, nspin, nlorentz
        complex(c_double_complex), dimension(len_ff), intent(in) :: ff
        type(fourvector) :: qvec
        type(fourvector), dimension(nin), intent(in) :: mom_in
        type(fourvector), dimension(nout), intent(in) :: mom_out
        integer(c_int), dimension(nin), intent(in) :: pids_in
        integer(c_int), dimension(nout), intent(in) :: pids_out
        complex(c_double_complex), dimension(nspin, nlorentz), intent(out) :: cur

        call current_init_had(mom_in,mom_out,qvec)
        call define_spinors()
        J_mu(nspin_f*nspin_i,4) ! that's how I want the output
        allocate(J_mu(nspin_f,nspin_in,4)) !...last two indices are for proton and neutron, see how to obtain particle spins and id. 
      !....
      ff1(1)=ff(1)
      ff2(1)=ff(2)
      ffa(1)=ff(3)
      ffa(2)=0.0d0
      !....
      call det_Ja(ff1,ff2,ffa)
      call hadr_curr_matrix_el(J_mu)



        cur(1, 1) = 1
        cur(1, 2) = 2
        cur(1, 3) = 3
        cur(1, 4) = 4
        cur(2, 1) = 5
        cur(2, 2) = 6
        cur(2, 3) = 7
        cur(2, 4) = 8
        cur(3, 1) = 9
        cur(3, 2) = 10
        cur(3, 3) = 11
        cur(3, 4) = 12
        cur(4, 1) = 13
        cur(4, 2) = 14
        cur(4, 3) = 15
        cur(4, 4) = 16
    end subroutine




   subroutine test_currents(self, mom_in, nin, mom_out, nout, qvec, ff, len_ff, cur, len_cur)
        use iso_c_binding
        use libvectors

        use dirac_matrices
        
        
        class(test), intent(inout) :: self
        integer(c_size_t), intent(in), value :: nin, nout, len_ff
        integer(c_int), intent(out) :: len_cur
        complex(c_double_complex), dimension(len_ff), intent(in) :: ff
        type(fourvector) :: qvec
        type(fourvector), dimension(nin), intent(in) :: mom_in
        type(fourvector), dimension(nout), intent(in) :: mom_out
        complex(c_double_complex), dimension(:), intent(out), pointer :: cur

        len_cur = 8
        allocate(cur(len_cur))





      call current_init_had(mom_in,mom_out,qvec,nspin_in,nspin_f)
      call define_spinors()
      allocate(J_mu(nspin_f,nspin_in,4,2)) !...last two indices are for proton and neutron, see how to obtain particle spins and id. 
      !....
      ff1(1)=ff(1)
      ff2(1)=ff(2)
      ffa(1)=ff(3)
      ffa(2)=0.0d0
      !....
      call det_Ja(ff1,ff2,ffa)
      call hadr_curr_matrix_el(J_mu)

      ! not sure why I have 8 currents


       ! cur(1) = 1
       ! cur(2) = 2
       ! cur(3) = 3
       ! cur(4) = 4
       ! cur(5) = 5
       ! cur(6) = 6
       ! cur(7) = 7
       ! cur(8) = 8
    end subroutine


    

        
        
end module
                
