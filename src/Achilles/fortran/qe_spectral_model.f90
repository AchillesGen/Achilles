module qe_spectral_model_test
    use iso_c_binding
    !use nuclear_model_interface
    use nuclear_model
    use libspectral_function
    implicit none
    private
    public :: qe_spec, build_qe_spec

    type(spectral_function) :: spectral_p, spectral_n

    type, extends(model) :: qe_spec
        contains
            procedure :: init => qe_spec_init
!            procedure :: fill_nucleus => qe_spec_fill
!            procedure :: nspins => qe_spec_spins
!            procedure :: states => qe_spec_states
            procedure :: currents => qe_spec_currents
            procedure :: ps_name => qe_spec_ps
            procedure :: mode => qe_spec_mode
            procedure :: cleanup => qe_spec_cleanup
    end type

contains

    function qe_spec_init(self, filename)
        use libutilities
        use dirac_matrices
        
        class(qe_spec), intent(inout) :: self
        integer :: ios, i
        character(len=*), intent(in) :: filename
        character(len=200) :: string
        integer, parameter :: read_unit = 99
        logical :: qe_spec_init
        character(len=:), allocatable :: trim_string 
        integer*8 :: length


       

        open(unit=read_unit, file=trim(filename), iostat=ios)
        if( ios /= 0 ) then
            qe_spec_init = .false.
            return
        endif

        read(read_unit, '(A)', iostat=ios) string
        trim_string = trim(string)
        length=len(trim_string)
        
        spectral_p = spectral_function(trim_string)
        write(6,*) 'trim_string',trim_string,len(trim_string)
        

        read(read_unit, '(A)', iostat=ios) string
        trim_string = trim(string)
        length=len(trim_string)
        
        spectral_n = spectral_function(trim_string)
        
        qe_spec_init = .true.
        write(6,*) 'prova', spectral_n%normalization()
        write(6,*) 'prova1', spectral_p%call(10.0d0,22.5d0)
 

        !....does it go here?

        call init(constants) ! load constants
        call dirac_matrices_in(constants%mqe) 

    end function

    function build_qe_spec() !..you register the different things that cAN BE CREATED
        class(model), pointer :: build_qe_spec
        allocate(qe_spec :: build_qe_spec)
    end function build_qe_spec

    function qe_spec_mode(self) !..interaction mode: QE, MEC, RES
        class(qe_spec), intent(inout) :: self
        integer :: qe_spec_mode
        qe_spec_mode = 1
    end function

    function qe_spec_ps(self) !...how to generate the nucler model phase space: HadronicMapper.hh
        class(qe_spec), intent(inout) :: self
        character(len=:), allocatable :: qe_spec_ps
        qe_spec_ps = "QESpectral"
    end function

    subroutine qe_spec_cleanup(self)
        class(qe_spec), intent(inout) :: self
    end subroutine




    subroutine qe_spec_currents(self, pids_in, mom_in, nin, pids_out, mom_out, nout, qvec, ff, len_ff, cur, nspin, nlorentz)
        use iso_c_binding
        use libvectors
        use dirac_matrices
        use libutilities
        
        
        class(qe_spec), intent(inout) :: self
        integer(c_size_t), intent(in), value :: nin, nout, len_ff, nspin, nlorentz
        complex(c_double_complex), dimension(len_ff), intent(in) :: ff 
        type(fourvector) :: qvec
        type(fourvector), dimension(nin), intent(in) :: mom_in
        type(fourvector), dimension(nout), intent(in) :: mom_out
        integer(c_long), dimension(nin), intent(in) :: pids_in
        integer(c_long), dimension(nout), intent(in) :: pids_out
        complex(c_double_complex), dimension(nspin, nlorentz), intent(out) :: cur

        integer(c_size_t) :: i,j        
        complex(c_double_complex), dimension(2) :: ffa
        complex(c_double_complex) :: ff1,ff2
        double precision, dimension(4) :: p4,pp4,q4
        complex(c_double_complex), dimension(2,2, nlorentz) :: J_mu
        double precision :: pmom, E, pke 


        p4=mom_in(1)%to_array()
        pp4=mom_out(1)%to_array()
        q4=qvec%to_array()
        
        call current_init_had(p4,pp4,q4) 
        call define_spinors()
        ff1=ff(1)
        ff2=ff(2)
        ffa(1)=ff(3)
        ffa(2)=0.0d0

 
        call det_Ja(ff1,ff2,ffa)
        call hadr_curr_matrix_el(J_mu)
        cur=(0.0d0,0.0d0)
        pmom=sqrt(sum(p4(2:4)**2))
        
        E=-p4(1)+constants%mqe
        pke=sqrt(spectral_p%call(pmom,E))
        !write(6,*) 'check', p4(1:4)

        do i=1,2
           do j=1,2
              cur(i+2*(j-1),:)= J_mu(j,i,:)*sqrt(spectral_p%call(pmom,E))
            enddo   
        enddo
     return
    end subroutine


end module qe_spectral_model_test

