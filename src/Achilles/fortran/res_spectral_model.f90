module res_spectral_model_test
    use iso_c_binding
    !use nuclear_model_interface
    use nuclear_model
    use libspectral_function   !...spectral_function_mod.f90
    implicit none
    private
    public :: res_spec, build_res_spec

    type(spectral_function) :: spectral_p, spectral_n

    type, extends(model) :: res_spec
        contains
            procedure :: init => res_spec_init
!            procedure :: fill_nucleus => qe_spec_fill
!            procedure :: nspins => qe_spec_spins
!            procedure :: states => qe_spec_states
            procedure :: currents => res_spec_currents
            procedure :: ps_name => res_spec_ps
            procedure :: mode => res_spec_mode
            procedure :: cleanup => res_spec_cleanup
    end type

contains

    function res_spec_init(self, filename)
        use libutilities
        use dirac_matrices_pi
        
        class(res_spec), intent(inout) :: self
        integer :: ios, i
        character(len=*), intent(in) :: filename
        character(len=200) :: string
        integer, parameter :: read_unit = 99
        logical :: res_spec_init
        character(len=:), allocatable :: trim_string 
        integer*8 :: length


       

        open(unit=read_unit, file=trim(filename), iostat=ios)
        if( ios /= 0 ) then
            res_spec_init = .false.
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
        
        res_spec_init = .true.
        !write(6,*) 'prova', spectral_n%normalization()
        !write(6,*) 'prova1', spectral_p%call(10.0d0,22.5d0)
 


        call init(constants) ! load constants

        call dirac_matrices_in(constants%mp,constants%mn,constants%meta,constants%mpi0,constants%mpip,constants%pi,constants%hbarc)

    end function

    function build_res_spec() !..you register the different things that cAN BE CREATED
        class(model), pointer :: build_res_spec
        allocate(res_spec :: build_res_spec)
    end function build_res_spec

    function res_spec_mode(self) !..interaction mode: QE, MEC, RES..See NuclearModel.cc
        class(res_spec), intent(inout) :: self
        integer :: res_spec_mode
        res_spec_mode = 4 !..see enum in NuclearModel.hh
    end function

    function res_spec_ps(self) !...how to generate the nucler model phase space: HadronicMapper.hh, HadronicMapper.cc
        class(res_spec), intent(inout) :: self
        character(len=:), allocatable :: res_spec_ps
        res_spec_ps = "QESpectral"
    end function

    subroutine res_spec_cleanup(self)
        class(res_spec), intent(inout) :: self
    end subroutine




    subroutine res_spec_currents(self, pids_in, mom_in, nin, pids_out, mom_out, nout, qvec, ff, len_ff, cur, nspin, nlorentz)
        use iso_c_binding
        use libvectors
        use dirac_matrices_pi
        use libutilities
        
        
        class(res_spec), intent(inout) :: self
        integer(c_size_t), intent(in), value :: nin, nout, len_ff, nspin, nlorentz
        complex(c_double_complex), dimension(len_ff), intent(in) :: ff 
        type(fourvector) :: qvec
        type(fourvector), dimension(nin), intent(in) :: mom_in
        type(fourvector), dimension(nout), intent(in) :: mom_out
        integer(c_long), dimension(nin), intent(in) :: pids_in
        integer(c_long), dimension(nout), intent(in) :: pids_out
        complex(c_double_complex), dimension(nspin, nlorentz), intent(out) :: cur

        integer(c_size_t) :: i,j,mode      
        complex(c_double_complex), dimension(2) :: ffa
        complex(c_double_complex) :: ff1,ff2
        double precision, dimension(4) :: p4,pp4,kpi4,q4
        complex(c_double_complex), dimension(2,2, nlorentz) :: J_mu
        double precision :: pmom, E, pke, e_couple, sw, alpha, mN, Vud
        complex(c_double_complex) ::  ii, coupling

        !Include coupling constant achilles
        coupling = ff(5)

        p4=mom_in(1)%to_array()
        pp4=mom_out(1)%to_array()
        kpi4=mom_out(2)%to_array()
        
        q4=qvec%to_array()

        !Set hard coded kinematics
        !p4 = (/772.871731d0, 13.8063591d0, -74.2304546d0, 83.6718288d0/)
        !pp4 = (/1052.14827d0, 28.5519856d0, -19.3847977d0, 474.841456d0/)
        !kpi4 = (/159.169311d0, 22.1525867d0, -63.2254495d0, 51.2628948d0/)
        !q4 = pp4 + kpi4 - p4         !This is the true q
        
        !write(6,*)'p = ', p4
        !write(6,*)'pf = ', pp4 
        !write(6,*)'kpi = ', kpi4
        !write(6,*)'q = ', q4

        mN = (constants%mn + constants%mp)/2.0d0
        pmom=sqrt(sum(p4(2:4)**2))
        E=-p4(1)+mN
        pke=sqrt(spectral_p%call(pmom,E))

        call current_init(p4,pp4,q4,kpi4) 

        call hadr_curr_matrix_el(pids_in(1),pids_out(1),pids_out(2),ff,len_ff,J_mu)

        cur=(0.0d0,0.0d0)

        do i=1,2
           do j=1,2
              cur(i+2*(j-1),:)= coupling*J_mu(j,i,:)*sqrt(spectral_p%call(pmom,E))
            enddo   
        enddo
     return
    end subroutine


end module res_spectral_model_test


