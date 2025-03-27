module intf_spectral_model
    use iso_c_binding
    use nuclear_model
    use libspectral_function
    use liblogging
    implicit none
    private
    public :: intf_spec, build_intf_spec

    type(spectral_function) :: spectral_p_MF, spectral_n_MF

    type, extends(model) :: intf_spec
        contains
            procedure :: init => intf_spec_init
            procedure :: currents => intf_spec_currents
            procedure, nopass :: model_name => intf_spec_name
            procedure :: ps_name => intf_spec_ps
            procedure :: mode => intf_spec_mode
            procedure :: init_wgt => intf_spec_init_wgt
            procedure :: cleanup => intf_spec_cleanup
            procedure, nopass :: inspirehep => intf_inspirehep
    end type

    integer :: compute_1body = 1

contains

    function intf_spec_init(self, filename, params)
        use libutilities
        use libsystem
        use dirac_matrices_intf
        use libmap

        class(intf_spec), intent(inout) :: self
        integer :: ios, i
        character(len=*), intent(in) :: filename
        character(len=:), allocatable :: filepath
        type(map), intent(in) :: params
        character(len=200) :: string
        integer, parameter :: read_unit = 99
        logical :: intf_spec_init
        character(len=:), allocatable :: trim_string 
        integer*8 :: length
        character(len=256) :: error_message

        filepath = find_file(filename, "Interference Model")
        call logger%debug("Interference Model: Loading param file "//trim(filepath))
        open(unit=read_unit, file=trim(filepath), iostat=ios, iomsg=error_message, status='old')
        if( ios /= 0 ) then
            intf_spec_init = .false.
            call logger%error("Interference Model: "//error_message)
            close(read_unit)
            return
        endif

        read(read_unit, '(A)', iostat=ios, iomsg=error_message) string
        if( ios /= 0 ) then
            intf_spec_init = .false.
            call logger%error("Interference Model: "//error_message)
            close(read_unit)
            return
        endif
        trim_string = trim(string)
        length=len(trim_string)
        call logger%debug("Interference Model: Using proton spectral function file "//trim_string)
        spectral_p_MF = spectral_function(trim_string)

        read(read_unit, '(A)', iostat=ios, iomsg=error_message) string
        if( ios /= 0 ) then
            intf_spec_init = .false.
            call logger%error("Interference Model: "//error_message)
            close(read_unit)
            return
        endif
        trim_string = trim(string)
        length=len(trim_string)
        call logger%debug("Interference Model: Using proton spectral function file "//trim_string)
        spectral_n_MF = spectral_function(trim_string)
        intf_spec_init = .true.

        if(params%empty().eqv. .true.) then
            print*,'This model requires a parameter file. Please check your run card for ModelParamsFile.'
            intf_spec_init = .false.
        endif

        call init(constants) ! load constants
        call dirac_matrices_in(constants%mdelta,constants%mqe,constants%mpip,constants%mrho,params%lookup("fpind"),params%lookup("fstar"),params%lookup("fpinn2"),params%lookup("ga"),params%lookup("lpi"),params%lookup("lpind")) 
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

    function intf_inspirehep() !...reference for the model
        character(len=:), allocatable :: intf_inspirehep
        intf_inspirehep = "Lovato:2023khk" ! TODO: Add inspirehep information
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
        complex(c_double_complex), dimension(2) :: ffa
        complex(c_double_complex), dimension(2,2, nlorentz) :: J_mu_pi_dir, J_mu_del_dir
        complex(c_double_complex), dimension(2,2, nlorentz) :: J_mu_pi_exc, J_mu_del_exc
        complex(c_double_complex), dimension(2,2, nlorentz) :: J_mu_1b, J_mu
        logical :: has_axial

        ! Differentiate EM from CC & NC
        if(ff%lookup("FMecA5").eq.(0.0d0,0.0d0) .and. pids_in(1).eq.pids_out(1)) then
            has_axial = .false.
        else
            has_axial = .true.
        endif

        p1_4=mom_in(1)%to_array()
        p2_4=mom_spect(1)%to_array()
        pp1_4=mom_out(1)%to_array()
        pp2_4=mom_spect(1)%to_array()
        q4=qvec%to_array()

        ffa(1)=ff%lookup("FA")
        ffa(2)=ff%lookup("FAP")

        call current_init(p1_4,p2_4,pp1_4,pp2_4,q4,pids_in(1),pids_out(1),pids_spect(1),has_axial)
        call define_spinors()

        J_mu = (0.0d0,0.0d0)
        J_mu_1b = (0.0d0,0.0d0)
        J_mu_del_dir = (0.0d0,0.0d0)
        J_mu_pi_dir = (0.0d0,0.0d0)
        J_mu_del_exc = (0.0d0,0.0d0)
        J_mu_pi_exc = (0.0d0,0.0d0)

        ! Compute 1 body current first
        if(compute_1body.eq.1) then
            call det_J1(ff%lookup("F1"),ff%lookup("F2"),ffa)
            call onebody_curr_matrix_el(J_mu_1b)
            J_mu = J_mu_1b
            !Switch flag so we compute 2 body next time
            compute_1body = 0
        else
            !Compute 2 body current
            !Direct current pieces first
            err = det_JaJb_JcJd(ff%lookup("FMecV3"),ff%lookup("FMecV4"),ff%lookup("FMecV5"),ff%lookup("FMecA5"),1)
            ! Avoid interpolating outside
            ! of delta potential range
            if (err.eq.1) then
                cur=(0.0d0,0.0d0)
                return
            endif
            call det_Jpi(ff%lookup("FPiEM"));
            call twobody_del_curr_matrix_el(J_mu_del_dir)
            call twobody_pi_curr_matrix_el(J_mu_pi_dir)

            !Now exchange currents
            err = det_JaJb_JcJd(ff%lookup("FMecV3"),ff%lookup("FMecV4"),ff%lookup("FMecV5"),ff%lookup("FMecA5"),2)
            ! Avoid interpolating outside
            ! of delta potential range
            if (err.eq.1) then
                cur=(0.0d0,0.0d0)
                return
            endif
            call det_Jpi(ff%lookup("FPiEM"));
            call twobody_del_curr_matrix_el(J_mu_del_exc)
            call twobody_pi_curr_matrix_el(J_mu_pi_exc)


            J_mu = (J_mu_pi_dir + J_mu_del_dir + J_mu_del_exc + J_mu_pi_exc)/(2.0d0*p2_4(1))

            !Switch flag so we compute 1 body next time
            compute_1body = 1
        endif

        do i=1,2
           do j=1,2
              cur(i+2*(j-1),:)=J_mu(j,i,:)
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
        real(c_double) :: wgt, pmom,pspec_mom,E, hole_wgt, spect_wgt

        p1_4=mom_in(1)%to_array()
        E=-p1_4(1)+constants%mqe
        pmom=sqrt(sum(p1_4(2:4)**2))

        p2_4=mom_spect(1)%to_array()
        pspec_mom=sqrt(sum(p2_4(2:4)**2)) 

        if (pids_in(1).eq.2212) then
            hole_wgt=spectral_p_MF%normalization()*spectral_p_MF%call(pmom,E)
        else
            hole_wgt=spectral_n_MF%normalization()*spectral_n_MF%call(pmom,E)
        endif

        if (pids_spect(1).eq.2212) then
            spect_wgt=spectral_p_MF%normalization()*spectral_p_MF%call(pspec_mom) 
        else
            spect_wgt=spectral_n_MF%normalization()*spectral_n_MF%call(pspec_mom) 
        endif

        wgt = hole_wgt * spect_wgt

    end function intf_spec_init_wgt


end module intf_spectral_model

