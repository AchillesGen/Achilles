module qe_spectral_model
    use iso_c_binding
    use nuclear_model
    use libspectral_function
    implicit none
    private
    public :: qe_spec, build_qe_spec

    type(spectral_function) :: spectral_p, spectral_n

    type, extends(model) :: qe_spec
        contains
            procedure :: init => qe_spec_init
            procedure :: currents => qe_spec_currents
            procedure, nopass :: model_name => qe_spec_name
            procedure :: ps_name => qe_spec_ps
            procedure :: mode => qe_spec_mode
            procedure :: init_wgt => qe_spec_init_wgt
            procedure :: cleanup => qe_spec_cleanup
            procedure, nopass :: inspirehep => qe_inspirehep
    end type

contains

    function qe_spec_init(self, filename, params)
        use libutilities
        use dirac_matrices
        use libmap
        
        class(qe_spec), intent(inout) :: self
        integer :: ios, i
        character(len=*), intent(in) :: filename
        type(map), intent(in) :: params
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

        read(read_unit, '(A)', iostat=ios) string
        trim_string = trim(string)
        length=len(trim_string)
        spectral_n = spectral_function(trim_string)
        qe_spec_init = .true.

        call init(constants) ! load constants
        call dirac_matrices_in(constants%mqe) 

        close(read_unit)
    end function

    function build_qe_spec() !..you register the different things that cAN BE CREATED
        class(model), pointer :: build_qe_spec
        allocate(qe_spec :: build_qe_spec)
    end function build_qe_spec

    function qe_spec_mode(self) !..interaction mode: QE, MEC, RES
        class(qe_spec), intent(inout) :: self
        integer :: qe_spec_mode
        qe_spec_mode = 2
    end function

    function qe_spec_name() !...name of the model
        character(len=:), allocatable :: qe_spec_name
        qe_spec_name = "QE_Spectral_Func"
    end function

    function qe_inspirehep() !...reference for the model
        character(len=:), allocatable :: qe_inspirehep
        qe_inspirehep = "Rocco:2018mwt" ! TODO: Add inspirehep information
    end function

    function qe_spec_ps(self) !...how to generate the nucler model phase space: HadronicMapper.hh
        class(qe_spec), intent(inout) :: self
        character(len=:), allocatable :: qe_spec_ps
        qe_spec_ps = "OneBodySpectral"
    end function

    subroutine qe_spec_cleanup(self)
        class(qe_spec), intent(inout) :: self
    end subroutine

    subroutine qe_spec_currents(self, pids_in, mom_in, nin, pids_out, mom_out, nout, pids_spect, mom_spect, nspect, qvec, ff, cur, nspin, nlorentz)
        use iso_c_binding
        use libvectors
        use dirac_matrices
        use libutilities
        use libmap
        implicit none
        
        class(qe_spec), intent(inout) :: self
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

        integer(c_size_t) :: i,j        
        complex(c_double_complex), dimension(2) :: ffa
        double precision, dimension(4) :: p4,pp4,q4
        complex(c_double_complex), dimension(2,2, nlorentz) :: J_mu

        p4=mom_in(1)%to_array()
        pp4=mom_out(1)%to_array()
        q4=qvec%to_array()
        
        call current_init_had(p4,pp4,q4) 
        call define_spinors()
        ffa(1)=ff%lookup("FA")
        ffa(2)=ff%lookup("FAP")
 
        call det_Ja(ff%lookup("F1"),ff%lookup("F2"),ffa)
        call hadr_curr_matrix_el(J_mu)
        cur=(0.0d0,0.0d0)

        do i=1,2
           do j=1,2
              cur(i+2*(j-1),:) = J_mu(j,i,:)
            enddo   
        enddo
     return
    end subroutine

    function qe_spec_init_wgt(self, pids_in, mom_in, nin, pids_spect, mom_spect, nspect, nproton, nneutron) result(wgt)
        use iso_c_binding
        use libvectors
        use libutilities

        class(qe_spec), intent(inout) :: self
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
        
    end function qe_spec_init_wgt
end module qe_spectral_model
