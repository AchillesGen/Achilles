module test_fortran_inteference
    use intf_spectral_model

    private
    public :: test_interference
    type(intf_spec) :: model

contains

    subroutine test_interference(test_suite)
        use unit_test
        use libmap
        type(test_suite_type) :: test_suite
        integer(4) :: i,j
        real(8),dimension(4,4) :: rmunu_np, rmunu_nn, rmunu_pn, rmunu_pp
        real(8),dimension(4,4) :: noemi_rmunu_np, noemi_rmunu_nn, noemi_rmunu_pn, noemi_rmunu_pp
        real(8),dimension(4,4) :: rmunu_wgt_np, noemi_rmunu_wgt_np
        character(len=100) :: file
        logical :: valid
        integer(8) :: pid_in, pid_out, pid_spect
        type(map) :: params

        file = 'data/info_intf_C12.data'
        call params%init()
        call params%insert('fpind', 0.54d0)
        call params%insert('fstar', 2.13d0)
        call params%insert('fpinn2', 0.081*4.0*acos(-1d0))
        call params%insert('ga', 1.26d0)
        call params%insert('lpi', 1300.0d0)
        call params%insert('lpind', 1150.0d0)

        call test_case_create('Fortran Spectral Inteference model', test_suite)
        valid = model%init(file, params)

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

        noemi_rmunu_np(1,1) = 2.7740890726777619d-4
        noemi_rmunu_np(1,2) = -3.7965422690241157d-4
        noemi_rmunu_np(1,3) = -2.3231176568845108d-4
        noemi_rmunu_np(1,4) = 4.1061096428384364d-4
        noemi_rmunu_np(2,1) = -4.3161000558052597d-4
        noemi_rmunu_np(2,2) = 2.5764714497521469d-3
        noemi_rmunu_np(2,3) = 2.0220656479648386d-4
        noemi_rmunu_np(2,4) = -2.4244417210859243d-3
        noemi_rmunu_np(3,1) = -7.1436783649723686d-4
        noemi_rmunu_np(3,2) = 1.6268240965361402d-4
        noemi_rmunu_np(3,3) = 4.3096945527350689d-3
        noemi_rmunu_np(3,4) = 1.0765531068087364d-4
        noemi_rmunu_np(4,1) = 4.0234418013681471d-4
        noemi_rmunu_np(4,2) = -2.4126274543634810d-3
        noemi_rmunu_np(4,3) = 2.2363340682766516d-4
        noemi_rmunu_np(4,4) = 2.3185632439470393d-3

        noemi_rmunu_nn(1,1) = -3.3020181556906658d-6
        noemi_rmunu_nn(1,2) = -1.3001929821670444d-4
        noemi_rmunu_nn(1,3) = -5.9958819055489485d-4
        noemi_rmunu_nn(1,4) = 4.3404339239395113d-5
        noemi_rmunu_nn(2,1) = -6.2171820196583785d-5
        noemi_rmunu_nn(2,2) = 4.0821888540744415d-4
        noemi_rmunu_nn(2,3) = -1.3414578647709242d-4
        noemi_rmunu_nn(2,4) = -4.0158278257070409d-4
        noemi_rmunu_nn(3,1) = 1.4288406444307087d-5
        noemi_rmunu_nn(3,2) = -6.8658892892674422d-5
        noemi_rmunu_nn(3,3) = -2.7998584954695824d-4
        noemi_rmunu_nn(3,4) = 3.3064525168051652d-5
        noemi_rmunu_nn(4,1) = 5.5817444445407680d-5
        noemi_rmunu_nn(4,2) = -4.1834036577303654d-4
        noemi_rmunu_nn(4,3) = -1.3053010694973955d-4
        noemi_rmunu_nn(4,4) = 3.7691312475286227d-4

        noemi_rmunu_pn(1,1) = -1.7885469087056062d-3
        noemi_rmunu_pn(1,2) = -7.8227095342742073d-4
        noemi_rmunu_pn(1,3) = -2.3389389889856866d-4
        noemi_rmunu_pn(1,4) = 2.5052082954706289d-5
        noemi_rmunu_pn(2,1) = -2.0521391860787834d-4
        noemi_rmunu_pn(2,2) = 3.8534930856720922d-3
        noemi_rmunu_pn(2,3) = 2.7741445279762321d-4
        noemi_rmunu_pn(2,4) = -3.4703733478158249d-3
        noemi_rmunu_pn(3,1) = 5.0733821810885205d-3
        noemi_rmunu_pn(3,2) = 8.5627806800367043d-4
        noemi_rmunu_pn(3,3) = 6.0672535084139809d-3
        noemi_rmunu_pn(3,4) = 1.7820579186363615d-3
        noemi_rmunu_pn(4,1) = 1.2206241730577179d-4
        noemi_rmunu_pn(4,2) = -3.6152094763973706d-3
        noemi_rmunu_pn(4,3) = 3.5427731528894448d-4
        noemi_rmunu_pn(4,4) = 3.3032579363516269d-3

        noemi_rmunu_pp(1,1) = 3.5130429995087597d-4
        noemi_rmunu_pp(1,2) = -1.5649428658579787d-4
        noemi_rmunu_pp(1,3) = -9.0439040078329200d-4
        noemi_rmunu_pp(1,4) = 1.5847629614406479d-4
        noemi_rmunu_pp(2,1) = -5.2010974451786391d-5
        noemi_rmunu_pp(2,2) = 6.0767008856932065d-4
        noemi_rmunu_pp(2,3) = -2.0035988588265645d-4
        noemi_rmunu_pp(2,4) = -5.8325143403928457d-4
        noemi_rmunu_pp(3,1) = 3.8351027153701806d-5
        noemi_rmunu_pp(3,2) = -9.9801838469260902d-5
        noemi_rmunu_pp(3,3) = -4.1488281068150548d-4
        noemi_rmunu_pp(3,4) = 5.3449626473519358d-5
        noemi_rmunu_pp(4,1) = 1.7778803543962618d-4
        noemi_rmunu_pp(4,2) = -6.0906130238444099d-4
        noemi_rmunu_pp(4,3) = -1.9776522199896026d-4
        noemi_rmunu_pp(4,4) = 5.8255693288547657d-4

        do i=1,4
            do j=1,4
                call assert_approximate(rmunu_np(i,j), noemi_rmunu_np(i,j), eps=1d-6, suite=test_suite)
                call assert_approximate(rmunu_nn(i,j), noemi_rmunu_nn(i,j), eps=1d-6, suite=test_suite)
                call assert_approximate(rmunu_pn(i,j), noemi_rmunu_pn(i,j), eps=1d-6, suite=test_suite)
                call assert_approximate(rmunu_pp(i,j), noemi_rmunu_pp(i,j), eps=1d-6, suite=test_suite)
                !print*,rmunu_pp(i,j)
            enddo
        enddo

	end subroutine

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
        call ff%insert('FA', (0.0d0,0.0d0))
        call ff%insert('FAP',(0.0d0,0.0d0))
        call ff%insert('FPiEM', (0.47263563216390736d0,0.0d0))
        call ff%insert('FMecV3', (1.2997163009564785d0,0.0d0))
        call ff%insert('FMecV4', (0.0d0,0.0d0))
        call ff%insert('FMecV5', (0.0,0.0d0))
        call ff%insert('FMecA5', (0.0,0.0d0))

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

        mom_in = fourvector(p1_4(1), p1_4(2), p1_4(3), p1_4(4))
        mom_out = fourvector(pp1_4(1), pp1_4(2), pp1_4(3), pp1_4(4))
        mom_spect = fourvector(p2_4(1), p2_4(2), p2_4(3), p2_4(4))
        qvec = fourvector(q4(1), q4(2), q4(3), q4(4))
        call model%currents(pids_in, mom_in, nin, pids_out, mom_out, nout, pids_spect, mom_spect, nspect, qvec, ff, cur_2b, nspin, nlorentz)


        ! cur_1b and cur_2b are (nspin,nlorentz) arrays
        ! want to compute J2mu^* J1nu + J2mu * J1nu^*
        do mu=1,4
            do nu=1,4
                do k=1,4
                    test_rmunu(mu,nu) = test_rmunu(mu,nu) + V*( cur_2b(k,mu)*conjg(cur_1b(k,nu)) + conjg(cur_2b(k,mu))*cur_1b(k,nu) )
                enddo
            enddo
        enddo

        ! Include flux factor and spin average for initial struck nucleon, 
        ! And phase space factor for final struck nucleon
        test_rmunu(:,:) = test_rmunu(:,:)/2.0d0/(2.0d0*E_OS)/(2.0d0*pp1_4(1))

	end function


end module test_fortran_inteference
