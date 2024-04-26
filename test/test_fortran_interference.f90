module test_fortran_inteference
    use intf_spectral_model

    private
    public :: test_interference
    type(intf_spec) :: model

contains

    subroutine test_interference(test_suite)
        use unit_test
        type(test_suite_type) :: test_suite
        integer(4) :: i,j
        real(8),dimension(4,4) :: rmunu_np, rmunu_nn, rmunu_pn, rmunu_pp
        real(8),dimension(4,4) :: noemi_rmunu_np, noemi_rmunu_nn, noemi_rmunu_pn, noemi_rmunu_pp
        real(8),dimension(4,4) :: rmunu_wgt_np, noemi_rmunu_wgt_np
        character(len=100) :: file
        logical :: valid
        integer(8) :: pid_in, pid_out, pid_spect

        file = 'data/intf_info.data'

        call test_case_create('Fortran Spectral Inteference model', test_suite)
        valid = model%init(file)

        pid_in = 2112 
        pid_out = 2112  
        pid_spect = 2212 
        rmunu_np = REAL(test_rmunu(pid_in,pid_out,pid_spect))

        pid_in = 2112 
        pid_out = 2112  
        pid_spect = 2112 
        rmunu_nn = REAL(test_rmunu(pid_in,pid_out,pid_spect))

        pid_in = 2212 
        pid_out = 2212  
        pid_spect = 2112 
        rmunu_pn = REAL(test_rmunu(pid_in,pid_out,pid_spect))

        pid_in = 2212 
        pid_out = 2212  
        pid_spect = 2212 
        rmunu_pp = REAL(test_rmunu(pid_in,pid_out,pid_spect))
        
        noemi_rmunu_np(1,1) = 2.25815838515586696d-4
        noemi_rmunu_np(1,2) = -6.87370002410698133d-4
        noemi_rmunu_np(1,3) = -5.99488419437026144d-4
        noemi_rmunu_np(1,4) = 6.22185108846922503d-4
        noemi_rmunu_np(2,1) = -4.43029375385752976d-4
        noemi_rmunu_np(2,2) = 1.85328194993925098d-3
        noemi_rmunu_np(2,3) = -5.08792832468896216d-4
        noemi_rmunu_np(2,4) = -1.86944960118112465d-3
        noemi_rmunu_np(3,1) = -7.60031481356077381d-4
        noemi_rmunu_np(3,2) = 1.11445617660166314d-4
        noemi_rmunu_np(3,3) = 2.67544322729919076d-3
        noemi_rmunu_np(3,4) = -5.70110813669765925d-5
        noemi_rmunu_np(4,1) = 2.24556015875737082d-4
        noemi_rmunu_np(4,2) = -2.29011535216880942d-3
        noemi_rmunu_np(4,3) = -5.07161148396928106d-4
        noemi_rmunu_np(4,4) = 2.05860507876153790d-3

        noemi_rmunu_nn(1,1) = 6.34934721365693256d-5
        noemi_rmunu_nn(1,2) = -4.32706464586625993d-4
        noemi_rmunu_nn(1,3) = -3.33308021256430143d-4
        noemi_rmunu_nn(1,4) = 3.68434157510645914d-4
        noemi_rmunu_nn(2,1) = -1.53784729737972676d-4
        noemi_rmunu_nn(2,2) = 2.37693265227252625d-4
        noemi_rmunu_nn(2,3) = -5.28922899693415934d-4
        noemi_rmunu_nn(2,4) = -3.29816787706393088d-4
        noemi_rmunu_nn(3,1) = -1.68661001040936804d-4
        noemi_rmunu_nn(3,2) = 1.14092705365056136d-4
        noemi_rmunu_nn(3,3) = 7.69631782745370152d-4
        noemi_rmunu_nn(3,4) = -7.13012253414570832d-5
        noemi_rmunu_nn(4,1) = -2.41046266030563968d-5
        noemi_rmunu_nn(4,2) = -7.58006488426452366d-4
        noemi_rmunu_nn(4,3) = -5.95616438115318205d-4
        noemi_rmunu_nn(4,4) =  5.95055766752572500d-4

        noemi_rmunu_pn(1,1) = -7.94952970657571926d-4
        noemi_rmunu_pn(1,2) = -1.12984978442762251d-3
        noemi_rmunu_pn(1,3) = -8.30193700581179019d-4
        noemi_rmunu_pn(1,4) = 6.21152801382668807d-4
        noemi_rmunu_pn(2,1) = 1.14234559375595772d-3
        noemi_rmunu_pn(2,2) = 2.92110864275982188d-3
        noemi_rmunu_pn(2,3) = -8.41978178583187580d-4
        noemi_rmunu_pn(2,4) = -2.28847671125400340d-3
        noemi_rmunu_pn(3,1) = 5.88881884116841149d-3
        noemi_rmunu_pn(3,2) = 8.69191335773667395d-4
        noemi_rmunu_pn(3,3) = 3.60660307327154324d-3
        noemi_rmunu_pn(3,4) = 1.77231258421293265d-3
        noemi_rmunu_pn(4,1) = 1.63514631230398180d-3
        noemi_rmunu_pn(4,2) = -3.25565878406849658d-3
        noemi_rmunu_pn(4,3) = -8.14939334984962171d-4
        noemi_rmunu_pn(4,4) = 3.38903566785703474d-3

        noemi_rmunu_pp(1,1) = 1.75027342139253794d-4
        noemi_rmunu_pp(1,2) = -6.31723616145496761d-4
        noemi_rmunu_pp(1,3) = -4.96916635775257335d-4
        noemi_rmunu_pp(1,4) = 5.66281239583195590d-4
        noemi_rmunu_pp(2,1) = 1.08693022365901129d-3
        noemi_rmunu_pp(2,2) = 4.83504090320769904d-4
        noemi_rmunu_pp(2,3) = -8.47654172398376593d-4
        noemi_rmunu_pp(2,4) = -1.39733140939000690d-4
        noemi_rmunu_pp(3,1) = 1.13425238829498675d-3
        noemi_rmunu_pp(3,2) = 3.07698902134448382d-4
        noemi_rmunu_pp(3,3) = 1.06908651972650588d-3
        noemi_rmunu_pp(3,4) = 2.60892112484071843d-4
        noemi_rmunu_pp(4,1) = 1.48692845967974738d-3
        noemi_rmunu_pp(4,2) = -9.67983235189612094d-4
        noemi_rmunu_pp(4,3) = -9.56686364586706414d-4
        noemi_rmunu_pp(4,4) = 1.28299611001954448d-3

        do i=1,4
            do j=1,4
                call assert_approximate(rmunu_np(i,j), noemi_rmunu_np(i,j), eps=1d-6, suite=test_suite)
                call assert_approximate(rmunu_nn(i,j), noemi_rmunu_nn(i,j), eps=1d-6, suite=test_suite)
                call assert_approximate(rmunu_pn(i,j), noemi_rmunu_pn(i,j), eps=1d-6, suite=test_suite)
                call assert_approximate(rmunu_pp(i,j), noemi_rmunu_pp(i,j), eps=1d-6, suite=test_suite)
            enddo
        enddo


        pid_in = 2112 
        pid_out = 2112  
        pid_spect = 2212 
        rmunu_wgt_np = REAL(test_rmunu_initw(pid_in,pid_out,pid_spect))

        noemi_rmunu_wgt_np(1,1) = 9.51224575707878716d-21
        noemi_rmunu_wgt_np(1,2) = -2.41248632849571453d-20
        noemi_rmunu_wgt_np(1,3) = 4.17210795707627211d-21
        noemi_rmunu_wgt_np(1,4) = 2.36386739042572073d-22
        noemi_rmunu_wgt_np(2,1) = 4.58060925308119883d-20
        noemi_rmunu_wgt_np(2,2) = 2.44833301279298096d-19
        noemi_rmunu_wgt_np(2,3) = -3.18927724010466878d-20
        noemi_rmunu_wgt_np(2,4) = 1.10490488903503561d-21
        noemi_rmunu_wgt_np(3,1) = 2.37041039205160801d-20
        noemi_rmunu_wgt_np(3,2) = -3.26948380764474816d-20
        noemi_rmunu_wgt_np(3,3) = 2.92776259895052242d-19
        noemi_rmunu_wgt_np(3,4) = 7.70564428141935593d-22
        noemi_rmunu_wgt_np(4,1) = -1.65007828995250356d-20
        noemi_rmunu_wgt_np(4,2) = 2.51364296744723778d-20
        noemi_rmunu_wgt_np(4,3) = 3.89856678742178511d-20
        noemi_rmunu_wgt_np(4,4) = -3.80347616354105316d-22

        do i=1,4
            do j=1,4
                call assert_approximate(rmunu_wgt_np(i,j), noemi_rmunu_wgt_np(i,j), eps=1d-3, suite=test_suite)
            enddo
        enddo

	end subroutine

    function test_rmunu_initw(pid_in, pid_out, pid_spect)
        use iso_c_binding
        use libmap
        use libvectors
        use libutilities
        implicit none 

        integer(c_size_t) :: nin, nout, nspect, nspin, nlorentz, mu, nu, k, nproton, nneutron
        type(fourvector) :: qvec
        double precision, dimension(4) :: p1_4,pp1_4,p2_4,pp2_4,q4,q4tilde
        double precision, parameter :: pi=acos(-1.0d0), kf=225.0d0, A=12.0d0, xmn=938.91875d0
        double precision :: rho,V,E_OS,wgt
        type(complex_map) :: ff
        integer(c_long) :: pid_in, pid_out, pid_spect
        integer(c_long), dimension(1) :: pids_in
        integer(c_long), dimension(1) :: pids_out
        integer(c_long), dimension(1) :: pids_spect
        integer(c_long) :: prot, neut
        type(fourvector), dimension(1) :: mom_in
        type(fourvector), dimension(1) :: mom_out
        type(fourvector), dimension(1) :: mom_spect
        complex(c_double_complex), dimension(4, 4) :: cur_1b,cur_2b
        complex(c_double_complex), dimension(4, 4) :: test_rmunu_initw
        character(len=100) :: key
        complex(c_double_complex) :: f1,f2
        test_rmunu_initw = (0.0d0,0.0d0)

        rho = (kf**3)/(1.5d0 * pi**2)
        V = rho/A

        pids_in = (/pid_in/)
        pids_out = (/pid_out/)
        pids_spect = (/pid_spect/)

        nin=1
        nout=1
        nspect=1
        nspin=4
        nlorentz=4
        nproton=6
        nneutron=6

        if(pid_in.eq.2112) then
            f1 = (-2.3609083681090711d-2,0.0d0)
            f2 = (-0.83569829475886337d0,0.0d0)
        else 
            f1 = (0.51525507443241114d0,0.0d0)
            f2 = (0.74314889999419664d0,0.0d0)
        endif

        call ff%init()
        call ff%insert('F1', f1)
        call ff%insert('F2', f2)
        call ff%insert('FA', (0d0,0d0))
        call ff%insert('FPiEM', (0.38597498154029247d0,0d0))
        call ff%insert('FMecV3', (1.0628969706765468d0,0d0))
        call ff%insert('FMecV4', (-0.57386454973548828d0,0d0))
        call ff%insert('FMecV5', (0.16571919637582450d0,0d0))
        call ff%insert('FMecA5', (0d0,0d0))

        !Hard coded kinematics
        p1_4 = (/976.96899597764241d0,3.7860376451409494d0,6.8794784523530295d0,-269.88578824231035d0/)
        E_OS = p1_4(1)
        p2_4 = (/953.198d0, 112.0d0, 55.0d0, 107.0d0/)
        pp1_4 = (/991.41875000000005d0, 3.7860376451409494d0, -6.8794784523530295d0, 318.24684985899410d0/)
        pp2_4 = p2_4 
        q4 =  (/70.0d0,0.0d0,0.0d0 ,588.13263810130445d0/)
        q4tilde =  (/14.449754022357638d0,0.0d0,0.0d0,588.13263810130445d0/)

        p1_4(1) = E_OS + q4tilde(1) - q4(1)


        mom_in = fourvector(p1_4(1), p1_4(2), p1_4(3), p1_4(4))
        mom_out = fourvector(pp1_4(1), pp1_4(2), pp1_4(3), pp1_4(4))
        mom_spect = fourvector(p2_4(1), p2_4(2), p2_4(3), p2_4(4))
        qvec = fourvector(q4(1), q4(2), q4(3), q4(4))


        call model%currents(pids_in, mom_in, nin, pids_out, mom_out, nout, pids_spect, mom_spect, nspect, qvec, ff, cur_1b, nspin, nlorentz)
        call model%currents(pids_in, mom_in, nin, pids_out, mom_out, nout, pids_spect, mom_spect, nspect, qvec, ff, cur_2b, nspin, nlorentz)


        ! cur_1b and cur_2b are (nspin,nlorentz) arrays
        ! want to compute J2mu^* J1nu + J2mu * J1nu^*
        do mu=1,4
            do nu=1,4
                do k=1,4
                    test_rmunu_initw(mu,nu) = test_rmunu_initw(mu,nu) + V*(cur_2b(k,mu)*cur_1b(k,nu) + conjg(cur_2b(k,mu)*cur_1b(k,nu)))
                enddo
            enddo
        enddo

        ! Include flux factor and spin average for initial struck nucleon, 
        ! And phase space factor for final struck nucleon
        test_rmunu_initw(:,:) = test_rmunu_initw(:,:)/2.0d0/(2.0d0*E_OS)/(2.0d0*pp1_4(1))

        wgt = model%init_wgt(pids_in, mom_in, nin, pids_spect, mom_spect, nspect, nproton, nneutron)

        test_rmunu_initw(:,:) = test_rmunu_initw*wgt

    end function

    function test_rmunu(pid_in, pid_out, pid_spect)
        use iso_c_binding
        use libmap
        use libvectors
        use libutilities
        implicit none 

        integer(c_size_t) :: nin, nout, nspect, nspin, nlorentz, mu, nu, k
        type(fourvector) :: qvec
        double precision, dimension(4) :: p1_4,pp1_4,p2_4,pp2_4,q4
        double precision, parameter :: pi=acos(-1.0d0), kf=225.0d0, A=12.0d0, xmn=938.91875d0
        double precision :: rho,V,E_OS
        type(complex_map) :: ff
        integer(c_long) :: pid_in, pid_out, pid_spect
        integer(c_long), dimension(1) :: pids_in
        integer(c_long), dimension(1) :: pids_out
        integer(c_long), dimension(1) :: pids_spect
        integer(c_long) :: prot, neut
        type(fourvector), dimension(1) :: mom_in
        type(fourvector), dimension(1) :: mom_out
        type(fourvector), dimension(1) :: mom_spect
        complex(c_double_complex), dimension(4, 4) :: cur_1b,cur_2b
        complex(c_double_complex), dimension(4, 4) :: test_rmunu
        character(len=100) :: key
        complex(c_double_complex) :: f1,f2
        test_rmunu = (0.0d0,0.0d0)

        rho = (kf**3)/(1.5d0 * pi**2)
        V = rho/A

        pids_in = (/pid_in/)
        pids_out = (/pid_out/)
        pids_spect = (/pid_spect/)

        nin=1
        nout=1
        nspect=1
        nspin=4
        nlorentz=4

        if(pid_in.eq.2112) then
            !print*,'Using neutron EM form factors'
            f1 = (-.017340697711910772d0,0.0d0)
            f2 = (-0.99011610632251468d0,0.0d0)
        else 
            !print*,'Using proton EM form factors'
            f1 = (0.59137717615671226d0,0.0d0)
            f2 = (0.89838354961982203d0,0.0d0)
        endif

        call ff%init()
        call ff%insert('F1', f1)
        call ff%insert('F2', f2)
        call ff%insert('FA', (0d0,0d0))
        call ff%insert('FPiEM', (0.47263563216390736d0,0d0))
        call ff%insert('FMecV3', (1.2997163009564785d0,0d0))
        call ff%insert('FMecV4', (-0.70172474887901359d0,0d0))
        call ff%insert('FMecV5', (0.21845330592141088d0,0d0))
        call ff%insert('FMecA5', (0d0,0d0))

        !Hard coded kinematics
        p1_4 = (/892.28663457549249d0,-77.036892302107034d0,-26.857649369116537d0,67.942565322702237d0/)
        p2_4 = (/1002.007566138d0, 96.6849111519814d0, -193.77595264909456d0, 274.8702302143906d0/)
        pp1_4 = (/1090.4949377893222d0, 281.3919769883637d0, -74.60793248850116d0, 470.7861894114025d0/)
        pp2_4 = p2_4 
        q4 =  (/198.20830321382988d0,358.42886929047074d0,-47.750283119384619d0 ,402.84362408870038d0/)

        E_OS = sqrt(xmn**2 + sum(p1_4(2:4)**2))

        mom_in = fourvector(p1_4(1), p1_4(2), p1_4(3), p1_4(4))
        mom_out = fourvector(pp1_4(1), pp1_4(2), pp1_4(3), pp1_4(4))
        mom_spect = fourvector(p2_4(1), p2_4(2), p2_4(3), p2_4(4))
        qvec = fourvector(q4(1), q4(2), q4(3), q4(4))


        call model%currents(pids_in, mom_in, nin, pids_out, mom_out, nout, pids_spect, mom_spect, nspect, qvec, ff, cur_1b, nspin, nlorentz)
        call model%currents(pids_in, mom_in, nin, pids_out, mom_out, nout, pids_spect, mom_spect, nspect, qvec, ff, cur_2b, nspin, nlorentz)


        ! cur_1b and cur_2b are (nspin,nlorentz) arrays
        ! want to compute J2mu^* J1nu + J2mu * J1nu^*
        do mu=1,4
            do nu=1,4
                do k=1,4
                    test_rmunu(mu,nu) = test_rmunu(mu,nu) + V*(cur_2b(k,mu)*cur_1b(k,nu) + conjg(cur_2b(k,mu)*cur_1b(k,nu)))
                enddo
            enddo
        enddo

        ! Include flux factor and spin average for initial struck nucleon, 
        ! And phase space factor for final struck nucleon
        test_rmunu(:,:) = test_rmunu(:,:)/2.0d0/(2.0d0*E_OS)/(2.0d0*pp1_4(1))

	end function


end module test_fortran_inteference