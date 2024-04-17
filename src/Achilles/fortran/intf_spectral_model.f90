module intf_spectral_model
    use iso_c_binding
    use nuclear_model
    use libspectral_function
    implicit none
    private
    public :: intf_spec, build_intf_spec

    type(spectral_function) :: spectral_p_MF, spectral_p_bkgd, &
        & spectral_n_MF, spectral_n_bkgd

    type, extends(model) :: intf_spec
        contains
            procedure :: init => intf_spec_init
            procedure :: currents => intf_spec_currents
            procedure, nopass :: model_name => intf_spec_name
            procedure :: ps_name => intf_spec_ps
            procedure :: mode => intf_spec_mode
            procedure :: init_wgt => intf_spec_init_wgt
            procedure :: cleanup => intf_spec_cleanup
    end type

    integer :: compute_1body = 1

contains

    function intf_spec_init(self, filename)
        use libutilities
        use dirac_matrices_intf
        
        class(intf_spec), intent(inout) :: self
        integer :: ios, i
        character(len=*), intent(in) :: filename
        character(len=200) :: string
        integer, parameter :: read_unit = 99
        logical :: intf_spec_init
        character(len=:), allocatable :: trim_string 
        integer*8 :: length

        open(unit=read_unit, file=trim(filename), iostat=ios)
        if( ios /= 0 ) then
            intf_spec_init = .false.
            return
        endif

        read(read_unit, '(A)', iostat=ios) string
        trim_string = trim(string)
        length=len(trim_string)
        spectral_p_MF = spectral_function(trim_string)

        read(read_unit, '(A)', iostat=ios) string
        trim_string = trim(string)
        length=len(trim_string)
        spectral_n_MF = spectral_function(trim_string)

        read(read_unit, '(A)', iostat=ios) string
        trim_string = trim(string)
        length=len(trim_string)
        spectral_p_bkgd = spectral_function(trim_string)

        read(read_unit, '(A)', iostat=ios) string
        trim_string = trim(string)
        length=len(trim_string)
        spectral_n_bkgd = spectral_function(trim_string)
        intf_spec_init = .true.

        call init(constants) ! load constants
        call dirac_matrices_in(1232.0d0,constants%mqe,constants%mpi0) 

        close(read_unit)
    end function

    function build_intf_spec() !..you register the different things that cAN BE CREATED
        class(model), pointer :: build_intf_spec
        allocate(intf_spec :: build_intf_spec)
    end function build_intf_spec

    function intf_spec_mode(self) !..interaction mode: QE, MEC, RES
        class(intf_spec), intent(inout) :: self
        integer :: intf_spec_mode
        intf_spec_mode = 7
    end function

    function intf_spec_name() !...name of the model
        character(len=:), allocatable :: intf_spec_name
        intf_spec_name = "Intf_Spectral_Func"
    end function

    function intf_spec_ps(self) !...how to generate the nucler model phase space: HadronicMapper.hh
        class(intf_spec), intent(inout) :: self
        character(len=:), allocatable :: intf_spec_ps
        intf_spec_ps = "IntfSpectral"
    end function

    subroutine intf_spec_cleanup(self)
        class(intf_spec), intent(inout) :: self
    end subroutine

    subroutine intf_spec_currents(self, pids_in, mom_in, nin, pids_out, mom_out, nout, pids_spect, mom_spect, nspect, qvec, ff, cur, nspin, nlorentz)
        use iso_c_binding
        use libvectors
        use dirac_matrices_intf
        use libutilities
        use libmap
        implicit none
        
        class(intf_spec), intent(inout) :: self
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

        integer*4 :: err 
        integer(c_size_t) :: i,j     
        double precision, dimension(4) :: p1_4,pp1_4,p2_4,pp2_4,q4
        double precision :: rho,A,V
        complex(c_double_complex), dimension(2,2, nlorentz) :: J_mu_pi, J_mu_del
        complex(c_double_complex), dimension(2,2, nlorentz) :: J_mu_1b, J_mu

        rho = (225.0d0**3)/(1.5d0 * (3.141d0**2))
        A = 12.0d0
        V = rho/A

        p1_4=mom_in(1)%to_array()
        p2_4=mom_spect(1)%to_array()
        pp1_4=mom_out(1)%to_array()
        pp2_4=mom_spect(1)%to_array()
        q4=qvec%to_array()

        call current_init(p1_4,p2_4,pp1_4,pp2_4,q4,2,pids_in(1),pids_spect(1))
        call define_spinors()
        call det_Ja(ff%lookup("F1"),ff%lookup("F2"),ff%lookup("FA"))

        J_mu = (0.0d0,0.0d0)
        J_mu_1b = (0.0d0,0.0d0)
        J_mu_del = (0.0d0,0.0d0)
        J_mu_pi = (0.0d0,0.0d0)

        if(compute_1body.eq.1) then
            call onebody_curr_matrix_el(J_mu_1b)
            J_mu = J_mu_1b
            compute_1body = 0
        else
            call det_Jpi()
            ! Avoid interpolating outside
            ! of delta potential range
            err = det_JaJb_JcJd(ff%lookup("FMecV3"),ff%lookup("FMecV4"),ff%lookup("FMecV5"),ff%lookup("FMecA5"))
            if (err.eq.1) then
                cur=(0.0d0,0.0d0)
                return
            endif
            call twobody_curr_matrix_J1Jdel_exc(J_mu_del)
            call twobody_curr_matrix_J1Jpi_exc(J_mu_pi)
            !call onebody_curr_matrix_el(J_mu_1b)
            J_mu = -sqrt(V)*(J_mu_pi + J_mu_del)/(2.0d0*p2_4(1))
            !print*, J_mu
            !J_mu = J_mu_1b
            compute_1body = 1
        endif

        do i=1,2
           do j=1,2
              cur(i+2*(j-1),:)= J_mu(j,i,:)
            enddo   
        enddo
     return
    end subroutine

    function intf_spec_init_wgt(self, pids_in, mom_in, nin, pids_spect, mom_spect, nspect, nproton, nneutron) result(wgt)
        use iso_c_binding
        use libvectors
        use libutilities

        class(intf_spec), intent(inout) :: self
        integer(c_long), dimension(nin), intent(in) :: pids_in
        type(fourvector), dimension(nin), intent(in) :: mom_in
        type(fourvector), dimension(nspect), intent(in) :: mom_spect
        integer(c_long), dimension(nspect), intent(in) :: pids_spect
        integer(c_size_t), intent(in), value :: nin, nspect, nproton, nneutron
        double precision, dimension(4) :: p1_4,p2_4
        real(c_double) :: wgt, pmom,pspec_mom,E

        p1_4=mom_in(1)%to_array()
        E=-p1_4(1)+constants%mqe
        pmom=sqrt(sum(p1_4(2:4)**2))

        p2_4=mom_spect(1)%to_array()
        pspec_mom=sqrt(sum(p2_4(2:4)**2)) 

        if (pids_in(1) == 2212) then
            wgt=spectral_p_MF%normalization()*spectral_p_MF%call(pmom,E)
        else
            wgt=spectral_n_MF%normalization()*spectral_n_MF%call(pmom,E)
        endif

        if (pids_spect(1) == 2212) then
            wgt=wgt*spectral_p_MF%normalization()*spectral_p_MF%call(pmom) 
        else
            wgt=wgt*spectral_n_MF%normalization()*spectral_n_MF%call(pmom) 
        endif
        
    end function intf_spec_init_wgt


end module intf_spectral_model

