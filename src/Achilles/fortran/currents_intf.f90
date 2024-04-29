module dirac_matrices_intf
    use libinterpolate
    implicit none
    integer*4, private, save :: i_fl,np_del
    integer*4, private, parameter :: nspin_in=2,nspin_f=2
    complex*16, private, parameter :: czero = (0.0d0,0.0d0)
    complex*16, private, parameter :: cone  = (1.0d0,0.0d0)
    complex*16, private, parameter :: ci    = (0.0d0,1.0d0)
    real*8, private, parameter :: pi=acos(-1.0d0) 
    complex*16, private, save :: pi_elec_ff
    real*8, parameter :: lsq=0.71*1.e6   
    real*8, private, parameter :: fgnd=5.0d0,fpind=0.54d0
    real*8, private, parameter :: fstar=2.13d0,fpinn2=0.081*4.0d0*pi!1.0094d0! 2.14/2.13 from JUAN, !=0.08*4.0d0*pi ARTURO
    real*8, private, parameter :: ga=1.26d0
    integer*4, private, save :: ax 
    real*8, private, parameter :: lpi=1300.0d0,lpind=1150.0d0
    real*8, private, save :: mqe, qval
    complex*16, private, save :: cv3, cv4, cv5, ca5
    complex*16, private, save :: sig(3,2,2),id(2,2),id4(4,4),up(2),down(2)
    complex*16, allocatable, private, save :: up1(:,:),up2(:,:),upp1(:,:),upp2(:,:), &
            &   ubarp1(:,:),ubarp2(:,:),ubarpp1(:,:),ubarpp2(:,:)
    complex*16, private, save :: t1(2),t1p(2),t2(2)
    complex*16, private, save :: gamma_mu(4,4,5),g_munu(4,4), sigma_munu(4,4,4,4)
    complex*16, private, save :: p1_sl(4,4),p2_sl(4,4),pp1_sl(4,4),pp2_sl(4,4), &
         &   k1_sl(4,4),k2_sl(4,4),q_sl(4,4), &
         &   Pi_k1(4,4),Pi_k2(4,4),Pi_k1e(4,4),Pi_k2e(4,4)
    real*8, private, save ::  p1(4),p2(4),pp1(4),pp2(4),q(4),k1(4),k2(4)
    complex*16, private, save :: J_a_mu(4,4,4),J_b_mu(4,4,4),J_c_mu(4,4,4),J_d_mu(4,4,4)
    complex*16, private, save :: Je_a_mu(4,4,4),Je_b_mu(4,4,4),Je_c_mu(4,4,4),Je_d_mu(4,4,4)
    complex*16, private, save :: J_pif(4,4,4),J_sea1(4,4,4),J_sea2(4,4,4),J_pl1(4,4,4),J_pl2(4,4,4)
    complex*16, private, save :: J_1(4,4,4)    
    real*8, private, save :: xmd,xmn,xmpi,xmrho,w
    type(interp1d), private, save:: interp
    real*8, private, allocatable :: pdel(:),pot_del(:)
contains

subroutine dirac_matrices_in(xmd_in,xmn_in,xmpi_in,xmrho_in)
    implicit none
    integer*4 :: i,j
    real*8 :: xmd_in,xmn_in,xmpi_in,xmrho_in
    xmd=xmd_in
    xmn=xmn_in
    xmpi=xmpi_in
    xmrho=xmrho_in

    print*,xmd, xmn, xmpi, xmrho

    allocate(up1(nspin_in,4),upp1(nspin_f,4), &
            &   ubarp1(nspin_in,4),ubarpp1(nspin_f,4), &
            &   up2(nspin_in,4),upp2(nspin_f,4), &
            &   ubarp2(nspin_in,4),ubarpp2(nspin_f,4))

    sig(:,:,:)=czero
    id(:,:)=czero
    id(1,1)=cone;id(2,2)=cone
    sig(1,1,2)=cone;sig(1,2,1)=cone
    sig(2,1,2)=-ci;sig(2,2,1)=ci
    sig(3,1,1)=cone;sig(3,2,2)=-cone
    gamma_mu=czero
    gamma_mu(1:2,1:2,1)=id;gamma_mu(3:4,3:4,1)=-id
    id4=czero    
    id4(1:2,1:2)=id;id4(3:4,3:4)=id
    do i=2,4
      gamma_mu(1:2,3:4,i)=sig(i-1,:,:)
      gamma_mu(3:4,1:2,i)=-sig(i-1,:,:)
    enddo
    gamma_mu(1:2,3:4,5)=id
    gamma_mu(3:4,1:2,5)=id
    g_munu=czero
    g_munu(1,1)=cone;g_munu(2,2)=-cone;g_munu(3,3)=-cone;g_munu(4,4)=-cone
    up(1)=cone;up(2)=czero
    down(1)=czero;down(2)=cone
    do i=1,4
       do j=1,4
          sigma_munu(:,:,i,j)=ci*0.5d0*(matmul(gamma_mu(:,:,i),gamma_mu(:,:,j)) &
               &     -matmul(gamma_mu(:,:,j),gamma_mu(:,:,i)))
       enddo
    enddo

    ! Read in delta potential and
    ! set up 1D interpolation
    open(10, file='data/rho_0p5.dat')
    read(10,*) np_del
    allocate(pdel(np_del),pot_del(np_del))
    do i=1,np_del
        read(10,*) pdel(i),pot_del(i)
    enddo 
    interp = interp1d(pdel, pot_del, np_del, 1)
end subroutine 

subroutine define_spinors()
    implicit none
    integer*4 :: i
    complex*16 :: sigp1(2,2),sigp2(2,2),sigpp1(2,2),sigpp2(2,2)
    real*8 :: cp1,cp2,cpp1,cpp2
    sigp1=czero
    sigp2=czero
    sigpp1=czero
    sigpp2=czero
    !.....initialize quadrispinors
    up1=czero
    up2=czero
    upp1=czero
    upp2=czero
!.......initialize normalization factors
    cp1=sqrt((p1(1)+xmn))!/(2.0d0*p1(1)))
    cp2=sqrt((p2(1)+xmn))!/(2.0d0*p2(1)))
    cpp1=sqrt((pp1(1)+xmn))!/(2.0d0*pp1(1)))
    cpp2=sqrt((pp2(1)+xmn))!/(2.0d0*pp2(1)))
!.....define sigma*p
    do i=1,3
      sigp1=sigp1+sig(i,:,:)*p1(i+1)
      sigp2=sigp2+sig(i,:,:)*p2(i+1)
      sigpp1=sigpp1+sig(i,:,:)*pp1(i+1)
      sigpp2=sigpp2+sig(i,:,:)*pp2(i+1)
    enddo
!.....build quadri-spinors    
    up1(1,1:2)=up(:)
    up1(1,3:4)=matmul(sigp1(:,:),up(:))/(p1(1)+xmn)
    up1(2,1:2)=down(:)
    up1(2,3:4)=matmul(sigp1(:,:),down(:))/(p1(1)+xmn)
    up1(:,:)=cp1*up1(:,:)
!
    up2(1,1:2)=up(:)
    up2(1,3:4)=matmul(sigp2(:,:),up(:))/(p2(1)+xmn)
    up2(2,1:2)=down(:)
    up2(2,3:4)=matmul(sigp2(:,:),down(:))/(p2(1)+xmn)
    up2(:,:)=cp2*up2(:,:)
!
    upp1(1,1:2)=up(:)
    upp1(1,3:4)=matmul(sigpp1(:,:),up(:))/(pp1(1)+xmn)
    upp1(2,1:2)=down(:)
    upp1(2,3:4)=matmul(sigpp1(:,:),down(:))/(pp1(1)+xmn)
    upp1(:,:)=cpp1*upp1(:,:)
!
    upp2(1,1:2)=up(:)
    upp2(1,3:4)=matmul(sigpp2(:,:),up(:))/(pp2(1)+xmn)
    upp2(2,1:2)=down(:)
    upp2(2,3:4)=matmul(sigpp2(:,:),down(:))/(pp2(1)+xmn)
    upp2(:,:)=cpp2*upp2(:,:)
!
    ubarp1(1,1:2)=up(:)
    ubarp1(1,3:4)=-matmul(up(:),sigp1(:,:))/(p1(1)+xmn)
    ubarp1(2,1:2)=down(:)
    ubarp1(2,3:4)=-matmul(down(:),sigp1(:,:))/(p1(1)+xmn)
    ubarp1(:,:)=cp1*ubarp1(:,:)
!
    ubarp2(1,1:2)=up(:)
    ubarp2(1,3:4)=-matmul(up(:),sigp2(:,:))/(p2(1)+xmn)
    ubarp2(2,1:2)=down(:)
    ubarp2(2,3:4)=-matmul(down(:),sigp2(:,:))/(p2(1)+xmn)
    ubarp2(:,:)=cp2*ubarp2(:,:)
!
    ubarpp1(1,1:2)=up(:)
    ubarpp1(1,3:4)=-matmul(up(:),sigpp1(:,:))/(pp1(1)+xmn)
    ubarpp1(2,1:2)=down(:)
    ubarpp1(2,3:4)=-matmul(down(:),sigpp1(:,:))/(pp1(1)+xmn)
    ubarpp1(:,:)=cpp1*ubarpp1(:,:)
!
    ubarpp2(1,1:2)=up(:)
    ubarpp2(1,3:4)=-matmul(up(:),sigpp2(:,:))/(pp2(1)+xmn)
    ubarpp2(2,1:2)=down(:)
    ubarpp2(2,3:4)=-matmul(down(:),sigpp2(:,:))/(pp2(1)+xmn)
    ubarpp2(:,:)=cpp2*ubarpp2(:,:)
    return
end subroutine

subroutine current_init(p1_in,p2_in,pp1_in,pp2_in,q_in,i_fl_in,nuc1_pid_in,nuc1_pid_out,nuc2_pid_in,has_axial_in)
    implicit none
    integer*4 :: i,i_fl_in
    integer*8 :: nuc1_pid_in,nuc1_pid_out,nuc2_pid_in
    real*8 ::  p1_in(4),p2_in(4),pp1_in(4),pp2_in(4),q_in(4)
    logical :: has_axial_in

    ! 1 if neutrinos
    ! 0 else
    if(has_axial_in .eqv. .false.) then
        ax = 0  
    else
        ax = 1
    endif

    i_fl=i_fl_in

    !struck nucleon
    if(nuc1_pid_in.eq.2212) then
        !print*,'proton ph'
        t1 = up 
    else
        !print*,'neutron ph'
        t1 = down
    endif

    !knocked out nucleon
    if(nuc1_pid_out.eq.2212) then
        !print*,'proton ph'
        t1p = up 
    else
        !print*,'neutron ph'
        t1p = down
    endif

    !spectator nucleon
    if(nuc2_pid_in.eq.2212) then
        !print*,'proton spectator'
        t2 = up
    else
        !print*,'neutron spectator'
        t2 = down
    endif

    p1=p1_in
    pp1=pp1_in
    q=q_in

    w=q(1)
    q(1)=w+p1(1)
    p1(1)=sqrt(p1(2)**2+p1(3)**2+p1(4)**2+xmn**2) 
    q(1)=q(1)-p1(1)

    p2=p2_in
    pp2=pp2_in

    k1=pp2-p1
    k2=pp1-p2

    p1_sl=czero
    p2_sl=czero
    pp1_sl=czero
    pp2_sl=czero
    k1_sl=czero
    k2_sl=czero
    q_sl=czero
       
    do i=1,4
       p1_sl=p1_sl+g_munu(i,i)*gamma_mu(:,:,i)*p1(i)  
       p2_sl=p2_sl+g_munu(i,i)*gamma_mu(:,:,i)*p2(i)
       pp1_sl=pp1_sl+g_munu(i,i)*gamma_mu(:,:,i)*pp1(i)  
       pp2_sl=pp2_sl+g_munu(i,i)*gamma_mu(:,:,i)*pp2(i)
       k1_sl=k1_sl+g_munu(i,i)*gamma_mu(:,:,i)*k1(i)  
       k2_sl=k2_sl+g_munu(i,i)*gamma_mu(:,:,i)*k2(i)
       q_sl=q_sl+g_munu(i,i)*gamma_mu(:,:,i)*q(i)  
    enddo
    if(i_fl.eq.1) then
       Pi_k1(:,:)=matmul(gamma_mu(:,:,5),k1_sl(:,:))/(k1(1)**2-sum(k1(2:4)**2)-xmpi**2)
       Pi_k2(:,:)=matmul(gamma_mu(:,:,5),k2_sl(:,:))/(k2(1)**2-sum(k2(2:4)**2)-xmpi**2)
    elseif(i_fl.eq.2) then
       Pi_k1e(:,:)=matmul(gamma_mu(:,:,5),k1_sl(:,:))/(k1(1)**2-sum(k1(2:4)**2)-xmpi**2)
       Pi_k2e(:,:)=matmul(gamma_mu(:,:,5),k2_sl(:,:))/(k2(1)**2-sum(k2(2:4)**2)-xmpi**2)
    endif

    return
end subroutine

subroutine det_J1(f1v,f2v,fa)
  implicit none
  integer*4 :: mu,nu
  complex*16 :: f1v,f2v,fa  
 ! print*,f1v,f2v,fa

  do mu=1,4
     J_1(:,:,mu)=czero
     do nu=1,4
        J_1(:,:,mu)=J_1(:,:,mu)+ci*f2v*sigma_munu(:,:,mu,nu)&
             &     *g_munu(nu,nu)*q(nu)/2.0d0/xmn
     enddo
     J_1(:,:,mu)=J_1(:,:,mu)+f1v*gamma_mu(:,:,mu)
 enddo
 
 do mu=1,4
     J_1(:,:,mu)=J_1(:,:,mu)+fa*matmul(gamma_mu(:,:,mu),gamma_mu(:,:,5))
 enddo  

 
  return
end subroutine det_J1

subroutine onebody_curr_matrix_el(J_mu)
   implicit none
   integer*4 :: i1,f1,i
   complex*16 :: J_mu(nspin_f,nspin_in,4)

   do i1=1,2
    do f1=1,2
       do i=1,4
          J_mu(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_1(:,:,i),up1(i1,:)))
        enddo
     enddo
   enddo
   
   return
end subroutine

function det_JaJb_JcJd(cv3_in, cv4_in, cv5_in, ca5_in) result(err)
    implicit none
    integer*4 :: i,j,mu, err
    complex*16 :: cv3_in, cv4_in, cv5_in, ca5_in
    real*8 :: pa(4),pb(4),pc(4),pd(4),width,fpik1,fpik2,fpindk2,fpindk1
    real*8 :: gadelta,gbdelta,gcdelta,gddelta,pa2,pb2,pc2,pd2
    real*8 :: pot_pa,pot_pb,pot_pc,pot_pd
    complex*16 :: pa_sl(4,4),pb_sl(4,4),pc_sl(4,4),pd_sl(4,4)
    complex*16 :: xmd_a,xmd_b,xmd_c,xmd_d
    complex*16 :: j_a_1(4,4,4),j_a_2(4,4,4,4),RSa(4,4,4,4),RSb(4,4,4,4),j_b_1(4,4,4,4),j_b_2(4,4,4)
    complex*16 :: j_c_1(4,4,4),j_c_2(4,4,4,4),RSc(4,4,4,4),RSd(4,4,4,4),j_d_1(4,4,4,4),j_d_2(4,4,4)
    complex*16 :: J_a(4,4,4),J_b(4,4,4),J_c(4,4,4),J_d(4,4,4)
  
    cv3 = cv3_in
    cv4 = cv4_in
    cv5 = cv5_in
    ca5 = ca5_in

    if(i_fl.eq.1)then
        pa(:)=p1(:)+q(:)
        pb(:)=pp1(:)-q(:)
        pc(:)=p2(:)+q(:)
        pd(:)=pp2(:)-q(:)
    elseif(i_fl.eq.2) then
        pa(:)=p1(:)+q(:)
        pb(:)=pp2(:)-q(:)
        pc(:)=p2(:)+q(:)
        pd(:)=pp1(:)-q(:)
    endif

    pa2=pa(1)**2-sum(pa(2:4)**2)
    pb2=pb(1)**2-sum(pb(2:4)**2)
    pc2=pc(1)**2-sum(pc(2:4)**2)
    pd2=pd(1)**2-sum(pd(2:4)**2)

    call delta_potential(pa2,pot_pa)
    call delta_potential(pb2,pot_pb)
    call delta_potential(pc2,pot_pc)
    call delta_potential(pd2,pot_pd)
    
    pa_sl=czero
    pb_sl=czero
    pc_sl=czero
    pd_sl=czero
    call delta_se(pa(1)**2-sum(pa(2:4)**2),gadelta,pot_pa)
    xmd_a=xmd-0.5d0*ci*gadelta
    call delta_se(pb(1)**2-sum(pb(2:4)**2),gbdelta,pot_pb)
    xmd_b=xmd-0.5d0*ci*gbdelta
    call delta_se(pc(1)**2-sum(pc(2:4)**2),gcdelta,pot_pc)
    xmd_c=xmd-0.5d0*ci*gcdelta
    call delta_se(pd(1)**2-sum(pd(2:4)**2),gddelta,pot_pd)
    xmd_d=xmd-0.5d0*ci*gddelta
    fpik1=(lpi**2-xmpi**2)/(lpi**2-k1(1)**2+sum(k1(2:4)**2))
    fpik2=(lpi**2-xmpi**2)/(lpi**2-k2(1)**2+sum(k2(2:4)**2))
    fpindk1=lpind**2/(lpind**2-k1(1)**2+sum(k1(2:4)**2))
    fpindk2=lpind**2/(lpind**2-k2(1)**2+sum(k2(2:4)**2))
    do i=1,4
       pa_sl=pa_sl+g_munu(i,i)*gamma_mu(:,:,i)*pa(i)
       pb_sl=pb_sl+g_munu(i,i)*gamma_mu(:,:,i)*pb(i)
       pc_sl=pc_sl+g_munu(i,i)*gamma_mu(:,:,i)*pc(i)
       pd_sl=pd_sl+g_munu(i,i)*gamma_mu(:,:,i)*pd(i)
    enddo
! costruisco i primi due termini della corrente a 2 corpi corrispondenti ai diagrammi a,b,c e d questo e' un passaggio intermedio,
! l'espressione finale di tali correnti e' data da j_a_mu, j_b_mu, j_c_mu, j_d_mu
    do i=1,4
      j_a_1(:,:,i)=k2(i)*id4(:,:)
      j_b_2(:,:,i)=k2(i)*id4(:,:)
      j_c_1(:,:,i)=k1(i)*id4(:,:)
      j_d_2(:,:,i)=k1(i)*id4(:,:)
      !!!...I AM USING THE FULL 
      do j=1,4
         RSa(:,:,i,j)=matmul(pa_sl(:,:)+xmd*id4(:,:),g_munu(i,j)*id4(:,:)-matmul(gamma_mu(:,:,i),gamma_mu(:,:,j))/3.0d0- &
    &    2.0d0*pa(i)*pa(j)/3.0d0/xmd**2*id4(:,:)-(gamma_mu(:,:,i)*pa(j)-gamma_mu(:,:,j)*pa(i))/3.0d0/xmd) &
    &    *(1.0d0/(pa(1)**2-sum(pa(2:4)**2)-xmd_a**2))
!   &     *(pa2-xmd**2)/((pa2-xmd**2)**2+xmd**2*ga**2)

         RSb(:,:,i,j)=matmul(pb_sl(:,:)+xmd*id4(:,:),g_munu(i,j)*id4(:,:)-matmul(gamma_mu(:,:,i),gamma_mu(:,:,j))/3.0d0- &
    &    2.0d0*pb(i)*pb(j)/3.0d0/xmd**2*id4(:,:)-(gamma_mu(:,:,i)*pb(j)-gamma_mu(:,:,j)*pb(i))/3.0d0/xmd) &
    &    *(1.0d0/(pb(1)**2-sum(pb(2:4)**2)-xmd_b**2))
!   &     *(pb2-xmd**2)/((pb2-xmd**2)**2+xmd**2*gb**2)
         RSc(:,:,i,j)=matmul(pc_sl(:,:)+xmd*id4(:,:),g_munu(i,j)*id4(:,:)-matmul(gamma_mu(:,:,i),gamma_mu(:,:,j))/3.0d0- &
    &    2.0d0*pc(i)*pc(j)/3.0d0/xmd**2*id4(:,:)-(gamma_mu(:,:,i)*pc(j)-gamma_mu(:,:,j)*pc(i))/3.0d0/xmd) &
    &    *(1.0d0/(pc(1)**2-sum(pc(2:4)**2)-xmd_c**2))
!   &     *(pc2-xmd**2)/((pc2-xmd**2)**2+xmd**2*gc**2)
         RSd(:,:,i,j)=matmul(pd_sl(:,:)+xmd*id4(:,:),g_munu(i,j)*id4(:,:)-matmul(gamma_mu(:,:,i),gamma_mu(:,:,j))/3.0d0- &
    &    2.0d0*pd(i)*pd(j)/3.0d0/xmd**2*id4(:,:)-(gamma_mu(:,:,i)*pd(j)-gamma_mu(:,:,j)*pd(i))/3.0d0/xmd) &
    &    *(1.0d0/(pd(1)**2-sum(pd(2:4)**2)-xmd_d**2))
!   &     *(pd2-xmd**2)/((pd2-xmd**2)**2+xmd**2*gd**2)
!         J_a_2(:,:,i,j)=cv3*matmul(g_munu(i,j)*q_sl(:,:)-q(i)*gamma_mu(:,:,j),gamma_mu(:,:,5))+ca5*xmn*g_munu(i,j)*id4(:,:)
!         J_b_1(:,:,i,j)=cv3*matmul(gamma_mu(:,:,5),g_munu(j,i)*q_sl(:,:)-q(j)*gamma_mu(:,:,i))+ca5*xmn*g_munu(j,i)*id4(:,:)
!         J_c_2(:,:,i,j)=cv3*matmul(g_munu(i,j)*q_sl(:,:)-q(i)*gamma_mu(:,:,j),gamma_mu(:,:,5))+ca5*xmn*g_munu(i,j)*id4(:,:)
!         J_d_1(:,:,i,j)=cv3*matmul(gamma_mu(:,:,5),g_munu(j,i)*q_sl(:,:)-q(j)*gamma_mu(:,:,i))+ca5*xmn*g_munu(j,i)*id4(:,:)

         J_a_2(:,:,i,j)=cv3*matmul(g_munu(i,j)*q_sl(:,:)-q(i)*gamma_mu(:,:,j),gamma_mu(:,:,5))+(cv4/xmn+cv5/xmn)*(g_munu(i,j)* &
     &        (q(1)*pa(1) -q(3)*pa(3))-q(i)*pa(j)) &
         &  + ca5*xmn*g_munu(i,j)*id4(:,:)

         J_b_1(:,:,i,j)=cv3*matmul(g_munu(i,j)*q_sl(:,:)-q(i)*gamma_mu(:,:,j),gamma_mu(:,:,5))+(cv4/xmn+cv5/xmn)*(g_munu(i,j)* &
     &        (q(1)*pb(1) -q(3)*pb(3))-q(i)*pb(j)) &
         &  + ca5*xmn*g_munu(i,j)*id4(:,:)

         J_c_2(:,:,i,j)=cv3*matmul(g_munu(i,j)*q_sl(:,:)-q(i)*gamma_mu(:,:,j),gamma_mu(:,:,5))+(cv4/xmn+cv5/xmn)*(g_munu(i,j)* &
     &        (q(1)*pc(1)-q(3)*pc(3))-q(i)*pc(j)) &
         &  + ca5*xmn*g_munu(i,j)*id4(:,:)

         J_d_1(:,:,i,j)=cv3*matmul(g_munu(i,j)*q_sl(:,:)-q(i)*gamma_mu(:,:,j),gamma_mu(:,:,5))+(cv4/xmn+cv5/xmn)*(g_munu(i,j)* &
     &        (q(1)*pd(1)-q(3)*pd(3))-q(i)*pd(j)) &
         &  + ca5*xmn*g_munu(i,j)*id4(:,:)
         

!         J_b_1(:,:,i,j)=cv3*matmul(gamma_mu(:,:,5),g_munu(j,i)*q_sl(:,:)-q(j)*gamma_mu(:,:,i))+ca5*xmn*g_munu(j,i)*id4(:,:)
!         J_c_2(:,:,i,j)=cv3*matmul(g_munu(i,j)*q_sl(:,:)-q(i)*gamma_mu(:,:,j),gamma_mu(:,:,5))+ca5*xmn*g_munu(i,j)*id4(:,:)
!         J_d_1(:,:,i,j)=cv3*matmul(gamma_mu(:,:,5),g_munu(j,i)*q_sl(:,:)-q(j)*gamma_mu(:,:,i))+ca5*xmn*g_munu(j,i)*id4(:,:)



      enddo
    enddo
! costruisco Jmua, Jmub
   do mu=1,4
      J_a(:,:,mu)=czero
      J_b(:,:,mu)=czero
      J_c(:,:,mu)=czero
      J_d(:,:,mu)=czero
      do i=1,4
         do j=1,4
            J_a(:,:,mu)=J_a(:,:,mu)+matmul(J_a_1(:,:,i)*g_munu(i,i),matmul(RSa(:,:,i,j),g_munu(j,j)*J_a_2(:,:,j,mu)))
            J_b(:,:,mu)=J_b(:,:,mu)+matmul(J_b_1(:,:,mu,i)*g_munu(i,i),matmul(RSb(:,:,i,j),g_munu(j,j)*J_b_2(:,:,j))) 
            J_c(:,:,mu)=J_c(:,:,mu)+matmul(J_c_1(:,:,i)*g_munu(i,i),matmul(RSc(:,:,i,j),g_munu(j,j)*J_c_2(:,:,j,mu)))
            J_d(:,:,mu)=J_d(:,:,mu)+matmul(J_d_1(:,:,mu,i)*g_munu(i,i),matmul(RSd(:,:,i,j),g_munu(j,j)*J_d_2(:,:,j))) 
         enddo
      enddo
   enddo
   if(i_fl.eq.1) then
      J_a_mu=J_a*fpik2*fpindk2*sqrt(fpinn2)*fstar/xmpi**2/xmn
      J_b_mu=J_b*fpik2*fpindk2*sqrt(fpinn2)*fstar/xmpi**2/xmn
      J_c_mu=J_c*fpik1*fpindk1*sqrt(fpinn2)*fstar/xmpi**2/xmn
      J_d_mu=J_d*fpik1*fpindk1*sqrt(fpinn2)*fstar/xmpi**2/xmn
   elseif(i_fl.eq.2)then
      Je_a_mu=J_a*fpik2*fpindk2*sqrt(fpinn2)*fstar/xmpi**2/xmn
      Je_b_mu=J_b*fpik2*fpindk2*sqrt(fpinn2)*fstar/xmpi**2/xmn
      Je_c_mu=J_c*fpik1*fpindk1*sqrt(fpinn2)*fstar/xmpi**2/xmn
      Je_d_mu=J_d*fpik1*fpindk1*sqrt(fpinn2)*fstar/xmpi**2/xmn
   endif

   err = 0

end function

subroutine det_Jpi(pi_elec_ff_in)
   implicit none
   integer*4 :: mu
   complex*16 :: pi_elec_ff_in
   real*8 :: fpik1,fpik2,frho1,frho2,fact

   pi_elec_ff = pi_elec_ff_in

   fpik1=(lpi**2-xmpi**2)/(lpi**2-k1(1)**2+sum(k1(2:4)**2))
   fpik2=(lpi**2-xmpi**2)/(lpi**2-k2(1)**2+sum(k2(2:4)**2))
   frho1=1.0d0/(1.0d0 - (k1(1)**2-sum(k1(2:4)**2))/xmrho**2)
   frho2=1.0d0/(1.0d0 - (k2(1)**2-sum(k2(2:4)**2))/xmrho**2)

   !...this factor is needed to fulfill current conservation, see A3 Dekker
   fact=(k1(1)**2-sum(k1(2:4)**2)-xmpi**2)*(k2(1)**2-sum(k2(2:4)**2)-xmpi**2) &
        & *(1.0d0/(k1(1)**2-sum(k1(2:4)**2)-xmpi**2)/(k2(1)**2-sum(k2(2:4)**2)-xmpi**2) &
        & - 1.0d0/(k1(1)**2-sum(k1(2:4)**2)-xmpi**2)/(lpi**2-k1(1)**2+sum(k1(2:4)**2)) &
        & - 1.0d0/(k2(1)**2-sum(k2(2:4)**2)-xmpi**2)/(lpi**2-k2(1)**2+sum(k2(2:4)**2)))
   do mu=1,4
      J_pif(:,:,mu)=pi_elec_ff*(k1(mu)-k2(mu))*Pi_k1e(:,:)!*fact
      J_sea1(:,:,mu)=-pi_elec_ff*matmul(gamma_mu(:,:,5),gamma_mu(:,:,mu))-ax*frho1/ga*gamma_mu(:,:,mu)!/fpik2**2
      J_sea2(:,:,mu)=pi_elec_ff*matmul(gamma_mu(:,:,5),gamma_mu(:,:,mu))+ax*frho2/ga*gamma_mu(:,:,mu)!/fpik1**2
      J_pl1(:,:,mu)=ax*frho1/ga*q(mu)*q_sl(:,:)/(q(1)**2-q(4)**2-xmpi**2)
      J_pl2(:,:,mu)=-ax*frho2/ga*q(mu)*q_sl(:,:)/(q(1)**2-q(4)**2-xmpi**2)
   enddo
  J_pif=J_pif*fpik1*fpik2*fpinn2/xmpi**2 
  J_sea1=J_sea1*fpik2**2*fpinn2/xmpi**2
  J_sea2=J_sea2*fpik1**2*fpinn2/xmpi**2
  J_pl1=J_pl1*fpinn2/xmpi**2*fpik2**2
  J_pl2=J_pl2*fpinn2/xmpi**2*fpik1**2
  return
end subroutine 

subroutine twobody_curr_matrix_J1Jdel_exc(J_mu)
   implicit none
   integer*4 :: i1,f1,i2,f2,i,j
   complex*16 :: J_12a(2,2,4),J_12b(2,2,4)
   complex*16 :: Je_12a(2,2,4),Je_12b(2,2,4)   
   complex*16 :: Je_12a_dag(2,2,4),Je_12b_dag(2,2,4)
   complex*16 :: Je_12c(2,2,4),Je_12d(2,2,4)   
   complex*16 :: Je_12c_dag(2,2,4),Je_12d_dag(2,2,4)
   complex*16 :: Je_2(2,2),Je_2_dag(2,2)
   complex*16 :: Je_1(2,2),Je_1_dag(2,2)   
   complex*16 :: j1ja(nspin_f,nspin_in,4),j1jb(nspin_f,nspin_in,4)
   complex*16 :: j1jc(nspin_f,nspin_in,4),j1jd(nspin_f,nspin_in,4)
   complex*16 :: J_mu(nspin_f,nspin_in,4)
   complex*16 :: cta,ctb,ctc,ctd ! isospin coefficients

   do i1=1,2
      do f1=1,2
         Je_2(f1,i1)=sum(ubarpp1(f1,:)*matmul(Pi_k2e(:,:),up2(i1,:)))
         Je_2_dag(f1,i1)=conjg(Je_2(f1,i1))
         Je_1(f1,i1)=sum(ubarpp2(f1,:)*matmul(Pi_k1e(:,:),up1(i1,:)))
         Je_1_dag(f1,i1)=conjg(Je_1(f1,i1))


         do i=1,4
            !J_12a(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_a_mu(:,:,i),up1(i1,:)))
            
            Je_12a(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(Je_a_mu(:,:,i),up1(i1,:)))
            Je_12a_dag(f1,i1,i)=conjg(Je_12a(f1,i1,i))
            !
            Je_12b(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(Je_b_mu(:,:,i),up1(i1,:)))
            Je_12b_dag(f1,i1,i)=conjg(Je_12b(f1,i1,i))


            Je_12c(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(Je_c_mu(:,:,i),up2(i1,:)))
            Je_12c_dag(f1,i1,i)=conjg(Je_12c(f1,i1,i))

            Je_12d(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(Je_d_mu(:,:,i),up2(i1,:)))
            Je_12d_dag(f1,i1,i)=conjg(Je_12d(f1,i1,i))


         enddo
      enddo
   enddo   
   j1ja=czero
   j1jb=czero
   j1jc=czero
   j1jd=czero
   
   do i=1,4
    do i1=1,2
        do f1=1,2
            do i2=1,2
              j1jb(f1,i1,i)=j1jb(f1,i1,i)+Je_12b_dag(i2,i1,i)*Je_2_dag(f1,i2)
              j1jc(f1,i1,i)=j1jc(f1,i1,i)+Je_12c_dag(f1,i2,i)*Je_1_dag(i2,i1)
              !We may need ja and jd for exclusive case, even though they cancel for inclusive
              j1ja(f1,i1,i)=j1ja(f1,i1,i)+Je_12a_dag(i2,i1,i)*Je_2_dag(f1,i2)
              j1jd(f1,i1,i)=j1jd(f1,i1,i)+Je_12d_dag(f1,i2,i)*Je_1_dag(i2,i1)
            enddo
        enddo
    enddo
   enddo

   if(ax.eq.0) then
    cta = IDeltaA(t1,t2,t2,t1)
    ctb = IDeltaB(t1,t2,t2,t1)
    ctc = IDeltaC(t1,t2,t2,t1)
    ctd = IDeltaD(t1,t2,t2,t1)
   else
    cta = IDeltaAdag(t2,t1p,t1,t2)
    ctb = IDeltaBdag(t2,t1p,t1,t2)
    ctc = IDeltaCdag(t2,t1p,t1,t2)
    ctd = IDeltaDdag(t2,t1p,t1,t2)
    !print*,'cta = ', cta 
    !print*,'ctb = ', ctb 
    !print*,'ctc = ', ctc 
    !print*,'ctd = ', ctd
   endif

   J_mu=(j1ja(:,:,:)*cta + j1jb(:,:,:)*ctb +j1jc(:,:,:)*ctc + j1jd(:,:,:)*ctd)/2.0d0

   return
end subroutine twobody_curr_matrix_J1Jdel_exc


subroutine twobody_curr_matrix_J1Jpi_exc(J_mu)
   implicit none
   integer*4 :: i1,f1,i2,f2,i,j
   complex*16 :: J_12a(2,2,4),J_12b(2,2,4)
   complex*16 :: Je_12a(2,2,4),Je_12b(2,2,4)   
   complex*16 :: Je_12a_dag(2,2,4),Je_12b_dag(2,2,4)
   complex*16 :: Je_2(2,2),Je_2_dag(2,2)
   complex*16 :: Je_1(2,2),Je_1_dag(2,2)
   complex*16 :: Je_s2(2,2,4),Je_s2_dag(2,2,4)
   complex*16 :: Je_f(2,2,4),Je_f_dag(2,2,4),Je_s1(2,2,4),Je_s1_dag(2,2,4)
   complex*16 :: Je_p1(2,2,4),Je_p1_dag(2,2,4),Je_p2(2,2,4),Je_p2_dag(2,2,4)
   complex*16 :: j1jf(nspin_f,nspin_in,4),j1js(nspin_f,nspin_in,4),j1jp(nspin_f,nspin_in,4)
   complex*16 :: J_mu(nspin_f,nspin_in,4)
   complex*16 :: ctf,cts,ctp  
   
   do i1=1,2
      do f1=1,2
         Je_1(f1,i1)=sum(ubarpp2(f1,:)*matmul(Pi_k1e(:,:),up1(i1,:)))
         Je_1_dag(f1,i1)=conjg(Je_1(f1,i1))
         Je_2(f1,i1)=sum(ubarpp1(f1,:)*matmul(Pi_k2e(:,:),up2(i1,:)))
         Je_2_dag(f1,i1)=conjg(Je_2(f1,i1))
         do i=1,4
            Je_f(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(J_pif(:,:,i),up1(i1,:)))
            Je_f_dag(f1,i1,i)=conjg(Je_f(f1,i1,i))
            Je_s1(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(J_sea1(:,:,i),up1(i1,:)))
            Je_s1_dag(f1,i1,i)=conjg(Je_s1(f1,i1,i))
            Je_s2(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_sea2(:,:,i),up2(i1,:)))
            Je_s2_dag(f1,i1,i)=conjg(Je_s2(f1,i1,i))
            Je_p1(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(J_pl1(:,:,i),up1(i1,:)))
            Je_p1_dag(f1,i1,i)=conjg(Je_p1(f1,i1,i))
            Je_p2(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_pl2(:,:,i),up2(i1,:)))
            Je_p2_dag(f1,i1,i)=conjg(Je_p2(f1,i1,i))
         enddo
      enddo
   enddo   

   j1jf=czero
   j1js=czero
   j1jp=czero

   do i=1,4
    do i1=1,2
        do f1=1,2
            do i2=1,2
                j1jf(f1,i1,i)=j1jf(f1,i1,i)+Je_f_dag(i2,i1,i)*Je_2_dag(f1,i2)
                j1js(f1,i1,i)=j1js(f1,i1,i)+Je_s1_dag(i2,i1,i)*Je_2_dag(f1,i2) &
                          &  +Je_1_dag(i2,i1)*Je_s2_dag(f1,i2,i)
                !j1jp(f1,i1,i)=j1jp(f1,i1,i)+Je_p1_dag(i2,i1,i)*Je_2_dag(f1,i2) &
                          !&  +Je_1_dag(i2,i1)*Je_p2_dag(f1,i2,i)
            enddo
        enddo
    enddo
   enddo


   if(ax.eq.0) then
       ctf = -Iv(t1,t2,t2,t1) !Want to compute matrix element of (Iv^{dag} = -Iv)
       cts = -Iv(t1,t2,t2,t1) !Want to compute matrix element of (Iv^{dag} = -Iv)
       ctp = -Iv(t1,t2,t2,t1)
   else
       ctf = -Ivminus(t2,t1p,t1,t2) 
       cts = -Ivminus(t2,t1p,t1,t2) 
       ctp = -Ivminus(t2,t1p,t1,t2) 
   endif

   J_mu=(ctf*j1jf(:,:,:) + cts*j1js(:,:,:) + ctp*j1jp(:,:,:))/2.0d0

   return
 end subroutine twobody_curr_matrix_J1Jpi_exc
  
subroutine delta_se(pd2,width,pot)
   implicit none
   real*8 :: pd2,width,kpi,ekpi,eknuc,r2a,rfa,pot
   width=0.0d0
   if (pd2.ge.(xmpi+xmn)**2)then
      kpi=dsqrt(1.0d0/4.0d0/pd2*(pd2-(xmn+xmpi)**2)*(pd2-(xmn-xmpi)**2))
      ekpi=dsqrt(kpi**2+xmpi**2)
      eknuc=dsqrt(kpi**2+xmn**2)
      r2a= (eknuc-ekpi)**2 -4.0d0*kpi**2
      rfa=sqrt(0.95d0*xmn**2/(0.95d0*xmn**2-r2a))
!.....definizione Gamma
      width=(4.0d0*fpind)**2/12.0d0/pi/xmpi**2*kpi**3/sqrt(pd2)*(xmn+eknuc)  !&
!   &  *rfa**2!*(lpind**2/(lpind**2-xmpi**2))**2
!      width=120.0d0

!      width=0.38/(3.0d0*xmpi**2)*kpi**3/sqrt(pd2)*(xmn+eknuc)  &
!   &  *rfa**2-pot*2.0d0!*(lpind**2/(lpind**2-xmpi**2))**2


      return
   endif
!.....Eq.3.26 Phys.Rev. C 49 2650

   return
end subroutine

subroutine delta_potential(p2,delta_pot)
    implicit none
    real*8 :: p2, delta_pot
    real*8 :: p_interp

    p_interp = sqrt(p2) 

    if(p_interp.gt.interp%max()) then
        p_interp = interp%max()
    elseif (p_interp.lt.interp%min()) then
        p_interp = interp%min()
    endif

    delta_pot = interp%call(p_interp)
    return
end subroutine delta_potential

function Iv(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: Iv

    Iv = ci*(me(1,it1,itp1)*me(2,it2,itp2) - me(2,it1,itp1)*me(1,it2,itp2))

    return
end function Iv


function IDeltaA(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaA, c

    c = me(3,it2,itp2)

    IDeltaA = (2.*c/3.) + (Iv(it1,it2,itp1,itp2)/3.)

    return
end function IDeltaA

function IDeltaB(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaB, c

    c = me(3,it2,itp2) 

    IDeltaB = (2.*c/3.) - (Iv(it1,it2,itp1,itp2)/3.) 

    return
end function IDeltaB

function IDeltaC(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaC, c

    c = me(3,it1,itp1) 

    IDeltaC = (2.*c/3.) - (Iv(it1,it2,itp1,itp2)/3.)

    return
end function IDeltaC

function IDeltaD(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaD, c

    c = me(3,it1,itp1) 

    IDeltaD = (2.*c/3.) + (Iv(it1,it2,itp1,itp2)/3.) 

    return
end function IDeltaD

function Ivminus(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: Ivminus

    Ivminus = ci*( ( mev(2,it1,itp1)*mev(3,it2,itp2) - mev(3,it1,itp1)*mev(2,it2,itp2) ) & 
    &    - ci*( mev(3,it1,itp1)*mev(1,it2,itp2) - mev(1,it1,itp1)*mev(3,it2,itp2) ) )

    return
end function Ivminus

function IDeltaAdag(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaAdag, c

    c = (mev(1,it2,itp2) - ci*mev(2,it2,itp2))*iden(it1,itp1)

    IDeltaAdag = (2.*c/3.) + (Ivminus(it1,it2,itp1,itp2)/3.)

    return
end function IDeltaAdag

function IDeltaBdag(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaBdag, c

    c = (mev(1,it2,itp2) - ci*mev(2,it2,itp2))*iden(it1,itp1)

    IDeltaBdag = (2.*c/3.) - (Ivminus(it1,it2,itp1,itp2)/3.) 

    return
end function IDeltaBdag

function IDeltaCdag(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaCdag, c

    c = (mev(1,it1,itp1) - ci*mev(2,it1,itp1))*iden(it2,itp2)

    IDeltaCdag = (2.*c/3.) - (Ivminus(it1,it2,itp1,itp2)/3.)

    return
end function IDeltaCdag

function IDeltaDdag(it1,it2,itp1,itp2)
    implicit none
    complex*16 :: it1(2),it2(2),itp1(2),itp2(2)
    complex*16 :: IDeltaDdag, c

    c = (mev(1,it1,itp1) - ci*mev(2,it1,itp1))*iden(it2,itp2)

    IDeltaDdag = (2.*c/3.) + (Ivminus(it1,it2,itp1,itp2)/3.) 

    return
end function IDeltaDdag


function me(i,it,itp)
    implicit none
    integer*4 :: i
    complex*16 :: me, it(2),itp(2)

    me = sum(it(:)*matmul(sig(i,:,:),itp))

    return
end function me

function mev(i,it,itp)
    implicit none
    integer*4 :: i
    complex*16 :: mev, it(2),itp(2)

    mev = sum(itp(:)*matmul(sig(i,:,:),it))

    return
end function mev

function iden(it,itp)
    implicit none 
    complex*16:: iden, it(2),itp(2)

    iden = sum(itp(:)*matmul(id,it))
    return
end function iden


end module dirac_matrices_intf





