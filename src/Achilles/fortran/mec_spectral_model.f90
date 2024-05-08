module mec_spectral_model
    use iso_c_binding
    use nuclear_model
    use libspectral_function
    implicit none
    private
    public :: mec_spec, build_mec_spec

    type(spectral_function) :: spectral_p, spectral_n

    type, extends(model) :: mec_spec
        contains
            procedure :: init => mec_spec_init
            procedure :: currents => mec_spec_currents
            procedure, nopass :: model_name => mec_spec_name
            procedure :: ps_name => mec_spec_ps
            procedure :: mode => mec_spec_mode
            procedure :: init_wgt => mec_spec_init_wgt
            procedure :: cleanup => mec_spec_cleanup
    end type

contains

    function mec_spec_init(self, filename, params)
        use libutilities
        use dirac_matrices_mec
        use libmap

        class(mec_spec), intent(inout) :: self
        integer :: ios, i
        character(len=*), intent(in) :: filename
        type(map), intent(in) :: params
        character(len=200) :: string
        integer, parameter :: read_unit = 99
        logical :: mec_spec_init
        character(len=:), allocatable :: trim_string 
        integer*8 :: length

        open(unit=read_unit, file=trim(filename), iostat=ios)
        if( ios /= 0 ) then
            mec_spec_init = .false.
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
        mec_spec_init = .true.

        if(params%empty().eqv. .true.) then
            print*,'This model requires a parameter file. Please check your run card for ModelParamsFile.'
            mec_spec_init = .false.
        endif

        call init(constants) ! load constants
        call dirac_matrices_in(constants%mdelta,constants%mqe,constants%mpip,constants%mrho,params%lookup("fpind"),params%lookup("fstar"),params%lookup("fpinn2"),params%lookup("ga"),params%lookup("lpi"),params%lookup("lpind")) 
        close(read_unit)
    end function

    function build_mec_spec() !..you register the different things that cAN BE CREATED
        class(model), pointer :: build_mec_spec
        allocate(mec_spec :: build_mec_spec)
    end function build_mec_spec

    function mec_spec_mode(self) !..interaction mode: QE, MEC, RES
        class(mec_spec), intent(inout) :: self
        integer :: mec_spec_mode
        mec_spec_mode = 3
    end function

    function mec_spec_name() !...name of the model
        character(len=:), allocatable :: mec_spec_name
        mec_spec_name = "Mec_Spectral_Func"
    end function

    function mec_spec_ps(self) !...how to generate the nucler model phase space: HadronicMapper.hh
        class(mec_spec), intent(inout) :: self
        character(len=:), allocatable :: mec_spec_ps
        mec_spec_ps = "MecSpectral"
    end function

    subroutine mec_spec_cleanup(self)
        class(mec_spec), intent(inout) :: self
    end subroutine

    subroutine mec_spec_currents(self, pids_in, mom_in, nin, pids_out, mom_out, nout, pids_spect, mom_spect, nspect, qvec, ff, cur, nspin, nlorentz)
        use iso_c_binding
        use libvectors
        use dirac_matrices_mec
        use libutilities
        use libmap
        implicit none
        
        class(mec_spec), intent(inout) :: self
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
        integer(c_size_t) :: i,j,k,l     
        double precision, dimension(4) :: p1_4,pp1_4,p2_4,pp2_4,q4
        complex(c_double_complex), dimension(2,2,2,2, nlorentz) :: J_mu_pi_dir, J_mu_del_dir
        complex(c_double_complex), dimension(2,2,2,2, nlorentz) :: J_mu_pi_exc, J_mu_del_exc
        complex(c_double_complex), dimension(2,2,2,2, nlorentz) :: J_mu
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

        !call current_init(p1_4,p2_4,pp1_4,pp2_4,q4,pids_in(1),pids_out(1),pids_spect(1),has_axial)
        !call define_spinors()

        J_mu = (0.0d0,0.0d0)
        J_mu_del_dir = (0.0d0,0.0d0)
        J_mu_pi_dir = (0.0d0,0.0d0)
        J_mu_del_exc = (0.0d0,0.0d0)
        J_mu_pi_exc = (0.0d0,0.0d0)


        !Compute 2 body current
        !Direct current pieces first
        

        !Now exchange currents

        !Flatten fortran array 
        do i=1,2
            do j=1,2
                do k=1,2
                    do l=1,2
                        cur(i+2*(j-1)+4*(k-1)+8*(l-1),:)=J_mu(l,k,j,i,:)
                    enddo
                enddo
            enddo   
        enddo

     return
    end subroutine

    function mec_spec_init_wgt(self, pids_in, mom_in, nin, pids_spect, mom_spect, nspect, nproton, nneutron) result(wgt)
        use iso_c_binding
        use libvectors
        use libutilities

        class(mec_spec), intent(inout) :: self
        integer(c_long), dimension(nin), intent(in) :: pids_in
        type(fourvector), dimension(nin), intent(in) :: mom_in
        type(fourvector), dimension(nspect), intent(in) :: mom_spect
        integer(c_long), dimension(nspect), intent(in) :: pids_spect
        integer(c_size_t), intent(in), value :: nin, nspect, nproton, nneutron
        double precision, dimension(4) :: p1_4,p2_4
        real(c_double) :: wgt,pmom1,pmom2,E1,E2,hole1_wgt,hole2_wgt

        p1_4=mom_in(1)%to_array()
        E1=-p1_4(1)+constants%mqe
        pmom1=sqrt(sum(p1_4(2:4)**2))

        p2_4=mom_in(2)%to_array()
        E2=-p2_4(1)+constants%mqe
        pmom2=sqrt(sum(p2_4(2:4)**2)) 

        if (pids_in(1).eq.2212) then
            hole1_wgt=spectral_p%normalization()*spectral_p%call(pmom1,E1)
        else
            hole1_wgt=spectral_n%normalization()*spectral_n%call(pmom1,E1)
        endif

        if (pids_spect(1).eq.2212) then
            hole2_wgt=spectral_p%normalization()*spectral_p%call(pmom2,E2) 
        else
            hole2_wgt=spectral_n%normalization()*spectral_n%call(pmom2,E2) 
        endif

        wgt = hole1_wgt * hole2_wgt

    end function mec_spec_init_wgt


end module mec_spectral_model

