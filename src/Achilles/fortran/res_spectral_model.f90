module res_spectral_model
    use iso_c_binding
    use nuclear_model
    use libspectral_function  
    implicit none
    private
    public :: res_spec, build_res_spec

    type(spectral_function) :: spectral_p, spectral_n

    type, extends(model) :: res_spec
        contains
            procedure :: init => res_spec_init
            procedure :: currents => res_spec_currents
            procedure, nopass :: model_name => res_spec_name
            procedure :: ps_name => res_spec_ps
            procedure :: mode => res_spec_mode
            procedure :: frame => res_spec_frame
            procedure :: init_wgt => res_spec_init_wgt
            procedure :: cleanup => res_spec_cleanup
    end type

contains

    function res_spec_init(self, filename, params)
        use libutilities
        use dirac_matrices_pi
        use libmap
        
        class(res_spec), intent(inout) :: self
        integer :: ios, i
        character(len=*), intent(in) :: filename
        type(map), intent(in) :: params
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

        read(read_unit, '(A)', iostat=ios) string
        trim_string = trim(string)
        length=len(trim_string)
        spectral_n = spectral_function(trim_string)
        
        res_spec_init = .true.
        call init(constants) ! load constants
        call dirac_matrices_in(constants%mp,constants%mn,constants%meta,constants%mpi0,constants%mpip,constants%pi,constants%hbarc)

        close(read_unit)
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

    function res_spec_frame(self) !..frame: LAB, QZ
        class(res_spec), intent(inout) :: self
        integer :: res_spec_frame
        res_spec_frame = 1 !..see enum in NuclearModel.hh
    end function

    function res_spec_name() !..name of the model
        character(len=:), allocatable :: res_spec_name
        res_spec_name = "RES_Spectral_Func"
    end function

    function res_spec_ps(self) !...how to generate the nucler model phase space: HadronicMapper.hh, HadronicMapper.cc
        class(res_spec), intent(inout) :: self
        character(len=:), allocatable :: res_spec_ps
        res_spec_ps = "OneBodySpectral"
    end function

    subroutine res_spec_cleanup(self)
        class(res_spec), intent(inout) :: self
    end subroutine

    subroutine res_spec_currents(self, pids_in, mom_in, nin, pids_out, mom_out, nout, pids_spect, mom_spect, nspect, qvec, ff, cur, nspin, nlorentz)
        use iso_c_binding
        use libvectors
        use dirac_matrices_pi
        use libutilities
        use libmap
        
        class(res_spec), intent(inout) :: self
        integer(c_size_t), intent(in), value :: nin, nout, nspect, nspin, nlorentz
        type(complex_map), intent(in) :: ff
        type(fourvector) :: qvec
        integer(c_long), dimension(nin), intent(in) :: pids_in
        integer(c_long), dimension(nout), intent(in) :: pids_out
        integer(c_long), dimension(nspect), intent(in) :: pids_spect
        type(fourvector), dimension(nin), intent(in) :: mom_in
        type(fourvector), dimension(nout), intent(in) :: mom_out
        type(fourvector), dimension(nspect), intent(in) :: mom_spect
        complex(c_double_complex), dimension(nlorentz, nspin), intent(out) :: cur

        integer(c_size_t) :: i,j,mode      
        double precision, dimension(4) :: p4,pp4,kpi4,q4
        complex(c_double_complex), dimension(2,2, nlorentz) :: J_mu
        complex(c_double_complex) ::  coupling
        logical :: has_axial

        !Include coupling constant achilles
        coupling = ff%lookup("FResV")

        p4=mom_in(1)%to_array()
        pp4=mom_out(1)%to_array()
        kpi4=mom_out(2)%to_array()

        q4=qvec%to_array()

        if(ff%lookup("FResA") /= 0.0d0) then
            has_axial = .true.
        else
            has_axial = .false.
        endif

        call current_init(p4,pp4,q4,kpi4) 
        call hadr_curr_matrix_el(pids_in(1),pids_out(1),pids_out(2),has_axial,J_mu)
        cur=(0.0d0,0.0d0)

        do i=1,2
           do j=1,2
              cur(i+2*(j-1),:)= coupling*J_mu(j,i,:)
            enddo   
        enddo
     return
    end subroutine

    function res_spec_init_wgt(self, pids_in, mom_in, nin, pids_spect, mom_spect, nspect, nproton, nneutron) result(wgt)
        use iso_c_binding
        use libvectors
        use libutilities

        class(res_spec), intent(inout) :: self
        integer(c_long), dimension(nin), intent(in) :: pids_in
        type(fourvector), dimension(nin), intent(in) :: mom_in
        type(fourvector), dimension(nspect), intent(in) :: mom_spect
        integer(c_long), dimension(nspect), intent(in) :: pids_spect
        integer(c_size_t), intent(in), value :: nin, nspect, nproton, nneutron
        double precision, dimension(4) :: p4
        real(c_double) :: wgt, pmom, E

        p4=mom_in(1)%to_array()
        E=-p4(1)+constants%mqe
        pmom=sqrt(sum(p4(2:4)**2))

        if (pids_in(1) == 2212) then
            wgt=nproton*spectral_p%call(pmom,E)
        else
            wgt=nneutron*spectral_n%call(pmom,E)
        endif
    end function res_spec_init_wgt
end module res_spectral_model
