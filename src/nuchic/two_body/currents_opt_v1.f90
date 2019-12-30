module dirac_matrices
    implicit none
    integer*4, private, save :: i_fl
    complex*16, private, parameter :: czero = (0.0d0,0.0d0)
    complex*16, private, parameter :: cone  = (1.0d0,0.0d0)
    complex*16, private, parameter :: ci    = (0.0d0,1.0d0)
    real*8, private, parameter :: pi=acos(-1.0d0)    
    real*8, private, parameter :: fgnd=5.0d0,fpind=0.54d0
    real*8, private, parameter :: fstar=2.13d0, xmrho=775.8d0,ga=1.26d0,fpinn2=1.0094d0! 2.14/2.13 from JUAN, !=0.08*4.0d0*pi ARTURO
    real*8, private, parameter :: lpi=1300.0d0,lpind=1150.0d0
    real*8, private, save :: mqe, qval
    complex*16, private, save :: sig(3,2,2),id(2,2),id4(4,4),up(2),down(2)
    complex*16, private, save :: up1(2,4),up2(2,4),upp1(2,4),upp2(2,4), &
            &   ubarp1(2,4),ubarp2(2,4),ubarpp1(2,4),ubarpp2(2,4)
    complex*16, private, save :: gamma_mu(4,4,5),g_munu(4,4)
    complex*16, private, save :: p1_sl(4,4),p2_sl(4,4),pp1_sl(4,4),pp2_sl(4,4), &
         &   k1_sl(4,4),k2_sl(4,4),q_sl(4,4), &
         &   Pi_k1(4,4),Pi_k2(4,4),Pi_k1e(4,4),Pi_k2e(4,4)
    real*8, private, save ::  p1(4),p2(4),pp1(4),pp2(4),q(4),k1(4),k2(4)
    complex*16, private, save :: J_a_mu(4,4,4),J_b_mu(4,4,4),J_c_mu(4,4,4),J_d_mu(4,4,4)
    complex*16, private, save :: Je_a_mu(4,4,4),Je_b_mu(4,4,4),Je_c_mu(4,4,4),Je_d_mu(4,4,4)
    complex*16, private, save :: J_pif(4,4,4),J_sea1(4,4,4),J_sea2(4,4,4),J_pl1(4,4,4),J_pl2(4,4,4)    
    real*8, private,save :: xmd,xmn,xmpi
contains

subroutine dirac_matrices_in(xmd_in,xmn_in,xmpi_in)
    implicit none
    integer*4 :: i
    real*8 :: xmd_in,xmn_in,xmpi_in
    xmd=xmd_in
    xmn=xmn_in
    xmpi=xmpi_in
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
    cp1=sqrt((p1(1)+xmn)/(2.0d0*p1(1)))
    cp2=sqrt((p2(1)+xmn)/(2.0d0*p2(1)))
    cpp1=sqrt((pp1(1)+xmn)/(2.0d0*pp1(1)))
    cpp2=sqrt((pp2(1)+xmn)/(2.0d0*pp2(1)))
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

subroutine current_init(p1_in,p2_in,pp1_in,pp2_in,q_in,k1_in,k2_in,i_fl_in)
    implicit none
    integer*4 :: i,i_fl_in
    real*8 ::  p1_in(4),p2_in(4),pp1_in(4),pp2_in(4),q_in(4),k1_in(4),k2_in(4)
    i_fl=i_fl_in
    p1=p1_in
    p2=p2_in
    pp1=pp1_in
    pp2=pp2_in
    k1=k1_in
    k2=k2_in
    q=q_in
!
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

subroutine det_Jpi(gep)
   implicit none
   integer*4 :: mu
   real*8 :: fpik1,fpik2,gep,frho1,frho2,fact
   fpik1=(lpi**2-xmpi**2)/(lpi**2-k1(1)**2+sum(k1(2:4)**2))
   fpik2=(lpi**2-xmpi**2)/(lpi**2-k2(1)**2+sum(k2(2:4)**2))
   frho1=0.0d0!1.0d0/(1.0d0-(k1(1)**2-sum(k1(2:4)**2))/xmrho**2)
   frho2=0.0d0!1.0d0/(1.0d0-(k2(1)**2-sum(k2(2:4)**2))/xmrho**2)
   !...this factor is needed to fulfill current conservation, see A3 Dekker
   fact=(k1(1)**2-sum(k1(2:4)**2)-xmpi**2)*(k2(1)**2-sum(k2(2:4)**2)-xmpi**2) &
        & *(1.0d0/(k1(1)**2-sum(k1(2:4)**2)-xmpi**2)/(k2(1)**2-sum(k2(2:4)**2)-xmpi**2) &
        & - 1.0d0/(k1(1)**2-sum(k1(2:4)**2)-xmpi**2)/(lpi**2-k1(1)**2+sum(k1(2:4)**2)) &
        & - 1.0d0/(k2(1)**2-sum(k2(2:4)**2)-xmpi**2)/(lpi**2-k2(1)**2+sum(k2(2:4)**2)))
   do mu=1,4
      J_pif(:,:,mu)=gep*(k1(mu)-k2(mu))*Pi_k1(:,:)!*fact
      J_sea1(:,:,mu)=-gep*matmul(gamma_mu(:,:,5),gamma_mu(:,:,mu))-frho1/ga*gamma_mu(:,:,mu)!/fpik2**2
      J_sea2(:,:,mu)=gep*matmul(gamma_mu(:,:,5),gamma_mu(:,:,mu))+frho2/ga*gamma_mu(:,:,mu)!/fpik1**2
      J_pl1(:,:,mu)=frho1/ga*q(mu)*q_sl(:,:)/(q(1)**2-q(4)**2-xmpi**2)
      J_pl2(:,:,mu)=-frho2/ga*q(mu)*q_sl(:,:)/(q(1)**2-q(4)**2-xmpi**2)
   enddo
  J_pif=J_pif*fpik1*fpik2*fpinn2/xmpi**2 
  J_sea1=J_sea1*fpik2**2*fpinn2/xmpi**2
  J_sea2=J_sea2*fpik1**2*fpinn2/xmpi**2
  J_pl1=J_pl1*fpinn2/xmpi**2*fpik2**2
  J_pl2=J_pl2*fpinn2/xmpi**2*fpik1**2
  return
end subroutine  


subroutine det_JaJb_JcJd(cv3,ca5,np_del,pdel,pot_del)
    use mathtool
    implicit none
    integer*4 :: i,j,mu,np_del
    real*8 :: pa(4),pb(4),pc(4),pd(4),width,fpik1,fpik2,fpindk2,fpindk1
    real*8 :: ga,gb,gc,gd,cv3,ca5,pa2,pb2,pc2,pd2,pdel(np_del),pot_del(np_del)
    real*8 :: pot_pa,pot_pb,pot_pc,pot_pd
    complex*16 :: pa_sl(4,4),pb_sl(4,4),pc_sl(4,4),pd_sl(4,4)
    complex*16 :: xmd_a,xmd_b,xmd_c,xmd_d
    complex*16 :: j_a_1(4,4,4),j_a_2(4,4,4,4),RSa(4,4,4,4),RSb(4,4,4,4),j_b_1(4,4,4,4),j_b_2(4,4,4)
    complex*16 :: j_c_1(4,4,4),j_c_2(4,4,4,4),RSc(4,4,4,4),RSd(4,4,4,4),j_d_1(4,4,4,4),j_d_2(4,4,4)
    complex*16 :: J_a(4,4,4),J_b(4,4,4),J_c(4,4,4),J_d(4,4,4)
  
    
    pa(:)=p1(:)+q(:)
    pb(:)=pp1(:)-q(:)
    pc(:)=p2(:)+q(:)
    pd(:)=pp2(:)-q(:)
    pa2=sum(pa(2:4)**2)
    pb2=sum(pb(2:4)**2)
    pc2=sum(pc(2:4)**2)
    pd2=sum(pd(2:4)**2)
    
    pot_pa=0.0d0
    pot_pb=0.0d0
    pot_pc=0.0d0
    pot_pd=0.0d0
    if(sqrt(pa2).lt.pdel(np_del)) call interpolint(pdel,pot_del,np_del,sqrt(pa2),pot_pa,1)
    if(sqrt(pb2).lt.pdel(np_del)) call interpolint(pdel,pot_del,np_del,sqrt(pb2),pot_pb,1)
    if(sqrt(pc2).lt.pdel(np_del)) call interpolint(pdel,pot_del,np_del,sqrt(pc2),pot_pc,1)
    if(sqrt(pd2).lt.pdel(np_del)) call interpolint(pdel,pot_del,np_del,sqrt(pd2),pot_pd,1)
    
    pa_sl=czero
    pb_sl=czero
    pc_sl=czero
    pd_sl=czero
    call delta_se(pa(1)**2-sum(pa(2:4)**2),ga,pot_pa)
    xmd_a=xmd-0.5d0*ci*ga
    call delta_se(pb(1)**2-sum(pb(2:4)**2),gb,pot_pb)
    xmd_b=xmd-0.5d0*ci*gb
    call delta_se(pc(1)**2-sum(pc(2:4)**2),gc,pot_pc)
    xmd_c=xmd-0.5d0*ci*gc
    call delta_se(pd(1)**2-sum(pd(2:4)**2),gd,pot_pd)
    xmd_d=xmd-0.5d0*ci*gd
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
!   &     *(pa(1)**2-sum(pa(2:4)**2)-xmd**2)/((pa(1)**2-sum(pa(2:4)**2)-xmd**2)**2+xmd**2*ga**2)

         RSb(:,:,i,j)=matmul(pb_sl(:,:)+xmd*id4(:,:),g_munu(i,j)*id4(:,:)-matmul(gamma_mu(:,:,i),gamma_mu(:,:,j))/3.0d0- &
    &    2.0d0*pb(i)*pb(j)/3.0d0/xmd**2*id4(:,:)-(gamma_mu(:,:,i)*pb(j)-gamma_mu(:,:,j)*pb(i))/3.0d0/xmd) &
    &    *(1.0d0/(pb(1)**2-sum(pb(2:4)**2)-xmd_b**2))
!   &     *(pb(1)**2-sum(pb(2:4)**2)-xmd**2)/((pb(1)**2-sum(pb(2:4)**2)-xmd**2)**2+xmd**2*gb**2)
         RSc(:,:,i,j)=matmul(pc_sl(:,:)+xmd*id4(:,:),g_munu(i,j)*id4(:,:)-matmul(gamma_mu(:,:,i),gamma_mu(:,:,j))/3.0d0- &
    &    2.0d0*pc(i)*pc(j)/3.0d0/xmd**2*id4(:,:)-(gamma_mu(:,:,i)*pc(j)-gamma_mu(:,:,j)*pc(i))/3.0d0/xmd) &
    &    *(1.0d0/(pc(1)**2-sum(pc(2:4)**2)-xmd_c**2))
!   &     *(pc(1)**2-sum(pc(2:4)**2)-xmd**2)/((pc(1)**2-sum(pc(2:4)**2)-xmd**2)**2+xmd**2*gc**2)
         RSd(:,:,i,j)=matmul(pd_sl(:,:)+xmd*id4(:,:),g_munu(i,j)*id4(:,:)-matmul(gamma_mu(:,:,i),gamma_mu(:,:,j))/3.0d0- &
    &    2.0d0*pd(i)*pd(j)/3.0d0/xmd**2*id4(:,:)-(gamma_mu(:,:,i)*pd(j)-gamma_mu(:,:,j)*pd(i))/3.0d0/xmd) &
    &    *(1.0d0/(pd(1)**2-sum(pd(2:4)**2)-xmd_d**2))
!   &     *(pd(1)**2-sum(pd(2:4)**2)-xmd**2)/((pd(1)**2-sum(pd(2:4)**2)-xmd**2)**2+xmd**2*gd**2)
         J_a_2(:,:,i,j)=cv3*matmul(g_munu(i,j)*q_sl(:,:)-q(i)*gamma_mu(:,:,j),gamma_mu(:,:,5))+ca5*xmn*g_munu(i,j)*id4(:,:)
         J_b_1(:,:,i,j)=cv3*matmul(gamma_mu(:,:,5),g_munu(j,i)*q_sl(:,:)-q(j)*gamma_mu(:,:,i))+ca5*xmn*g_munu(j,i)*id4(:,:)
         J_c_2(:,:,i,j)=cv3*matmul(g_munu(i,j)*q_sl(:,:)-q(i)*gamma_mu(:,:,j),gamma_mu(:,:,5))+ca5*xmn*g_munu(i,j)*id4(:,:)
         J_d_1(:,:,i,j)=cv3*matmul(gamma_mu(:,:,5),g_munu(j,i)*q_sl(:,:)-q(j)*gamma_mu(:,:,i))+ca5*xmn*g_munu(j,i)*id4(:,:)
!         J_a_2(:,:,i,j)=0.5d0*cv3*matmul(q(j)*gamma_mu(:,:,i)-matmul(gamma_mu(:,:,j),matmul(q_sl(:,:),gamma_mu(:,:,i))), &
!    &                   gamma_mu(:,:,5))+ca5*xmn*g_munu(i,j)*id4(:,:)
!         J_b_1(:,:,i,j)=0.5d0*cv3*matmul(gamma_mu(:,:,5),q(i)*gamma_mu(:,:,j)-matmul(gamma_mu(:,:,j),matmul(q_sl(:,:),& 
!    &                   gamma_mu(:,:,i))))+ca5*xmn*g_munu(j,i)*id4(:,:)
!         J_c_2(:,:,i,j)=0.5d0*cv3*matmul(q(j)*gamma_mu(:,:,i)-matmul(gamma_mu(:,:,j),matmul(q_sl(:,:),gamma_mu(:,:,i))), &
!    &                   gamma_mu(:,:,5))+ca5*xmn*g_munu(i,j)*id4(:,:)
!         J_d_1(:,:,i,j)=0.5d0*cv3*matmul(gamma_mu(:,:,5),q(i)*gamma_mu(:,:,j)-matmul(gamma_mu(:,:,j),matmul(q_sl(:,:),& 
!    &                   gamma_mu(:,:,i))))+ca5*xmn*g_munu(j,i)*id4(:,:)


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

end subroutine


subroutine det_JaJc_dir(r_cc,r_cl,r_ll,r_t,r_tp)
   implicit none
   integer*4 :: i1,f1,i,j
   complex*16 ::tr_jaja(4,4),tr_2,tr_jbjb(4,4),tr_jajb(4,4)
   complex*16 ::tr_jc(4),tr_jd(4),tr_ja(4),tr_jb(4)
   complex*16 :: J_12a(2,2,4),J_12b(2,2,4),J_2(2,2)
   complex*16 :: J_12c(2,2,4),J_12d(2,2,4),J_1(2,2)
   complex*16 :: J_12a_dag(2,2,4),J_12b_dag(2,2,4),J_2_dag(2,2)
   complex*16 :: J_12c_dag(2,2,4),J_12d_dag(2,2,4),J_1_dag(2,2)
   complex*16 :: jaja(4,4),jbjb(4,4),jajb(4,4)
   complex*16 :: jajc(4,4),jajd(4,4),jbjc(4,4),jbjd(4,4),r_av(4,4)
   real*8 :: r_cc,r_cl,r_ll,r_t,r_tp
   


   J_12a=czero
   J_12b=czero
   J_12a_dag=czero
   J_12b_dag=czero
   J_2=czero
   J_2_dag=czero   
   J_12c=czero
   J_12d=czero
   J_12c_dag=czero
   J_12d_dag=czero
   J_1=czero
   J_1_dag=czero   

   
   tr_2=0.0d0
   tr_jc=0.0d0
   tr_jd=0.0d0
   tr_jaja=0.0d0
   tr_jbjb=0.0d0
   tr_jajb=0.0d0
   tr_ja=0.0d0
   tr_jb=0.0d0

   do i1=1,2
      do f1=1,2
         J_1(f1,i1)=sum(ubarpp1(f1,:)*matmul(Pi_k1(:,:),up1(i1,:)))
         J_1_dag(f1,i1)=conjg(J_1(f1,i1))
         J_2(f1,i1)=sum(ubarpp2(f1,:)*matmul(Pi_k2(:,:),up2(i1,:)))
         J_2_dag(f1,i1)=conjg(J_2(f1,i1))
         tr_2=tr_2+J_2_dag(f1,i1)*J_2(f1,i1)
         do i=1,4
           J_12a(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_a_mu(:,:,i),up1(i1,:)))
           J_12b(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_b_mu(:,:,i),up1(i1,:)))
           J_12a_dag(f1,i1,i)=conjg(J_12a(f1,i1,i))
           J_12b_dag(f1,i1,i)=conjg(J_12b(f1,i1,i))
           !!
           J_12c(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(J_c_mu(:,:,i),up2(i1,:)))
           J_12d(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(J_d_mu(:,:,i),up2(i1,:)))
         enddo
      enddo
   enddo   

   do i1=1,2
      do f1=1,2
         do i=1,4
            tr_ja(i)=tr_ja(i)+J_12a_dag(f1,i1,i)*J_1(f1,i1) !particle 1
            tr_jb(i)=tr_jb(i)+J_12b_dag(f1,i1,i)*J_1(f1,i1) !particle 1
            tr_jc(i)=tr_jc(i)+J_2_dag(f1,i1)*J_12c(f1,i1,i) !particle 2
            tr_jd(i)=tr_jd(i)+J_2_dag(f1,i1)*J_12d(f1,i1,i) !particle 2
            do j=1,4
               tr_jaja(i,j)=tr_jaja(i,j)+J_12a_dag(f1,i1,i)*J_12a(f1,i1,j)
               tr_jbjb(i,j)=tr_jbjb(i,j)+J_12b_dag(f1,i1,i)*J_12b(f1,i1,j)
               tr_jajb(i,j)=tr_jajb(i,j)+J_12a_dag(f1,i1,i)*J_12b(f1,i1,j) &
                    & +J_12b_dag(f1,i1,i)*J_12a(f1,i1,j)
            enddo
         enddo
      enddo  
   enddo     

   jaja(:,:)=tr_jaja(:,:)*tr_2*16.0d0/3.0d0
   jbjb(:,:)=tr_jbjb(:,:)*tr_2*16.0d0/3.0d0
   jajb(:,:)=tr_jajb(:,:)*tr_2*(16.0d0/9.0d0)
   do i=1,4
      do j=1,4
      jajc(i,j)=tr_ja(i)*tr_jc(j)*(-16.0d0/9.0d0)
      jajd(i,j)=tr_ja(i)*tr_jd(j)*(16.0d0/9.0d0)
      jbjc(i,j)=tr_jb(i)*tr_jc(j)*(16.0d0/9.0d0)
      jbjd(i,j)=tr_jb(i)*tr_jd(j)*(-16.0d0/9.0d0)
      enddo
   enddo  
   r_av=jaja+jbjb+jajb+jajc+jajd+jbjc+jbjd
   r_cc=r_av(1,1)
   r_cl=-0.5d0*(r_av(1,4)+r_av(4,1))
   r_ll=r_av(4,4)
   r_t=r_av(2,2)+r_av(3,3)
   r_tp=-0.5d0*ci*(r_av(2,3)-r_av(3,2)) 
 
  return
end subroutine det_JaJc_dir





subroutine det_JaJc_exc(r_cc,r_cl,r_ll,r_t,r_tp)
   implicit none
   integer*4 :: i1,f1,i2,f2,i,j,fe
   complex*16 ::tr_jaja(4,4),tr_jbjb(4,4),tr_jajb(4,4)
   complex*16 ::tr_jajc(4,4),tr_jajd(4,4),tr_jbjc(4,4),tr_jbjd(4,4)
   complex*16 :: J_12a(2,2,4),J_12b(2,2,4),J_2(2,2)
   complex*16 :: J_12c(2,2,4),J_12d(2,2,4),J_1(2,2)
   complex*16 :: Je_12a(2,2,4),Je_12b(2,2,4),Je_2(2,2)   
   complex*16 :: Je_12a_dag(2,2,4),Je_12b_dag(2,2,4),Je_2_dag(2,2)
   complex*16 :: Je_1(2,2),Je_1_dag(2,2)
   complex*16 :: jaja(4,4),jbjb(4,4),jajb(4,4)
   complex*16 :: jajc(4,4),jajd(4,4),jbjc(4,4),jbjd(4,4),r_exc(4,4)
   real*8 :: r_cc,r_cl,r_ll,r_t,r_tp
   


   J_12a=czero
   J_12b=czero
   Je_12a_dag=czero
   Je_12b_dag=czero
   J_2=czero
   Je_2_dag=czero   
   J_12c=czero
   J_12d=czero
   J_1=czero
   Je_1_dag=czero   

   
   tr_jaja=0.0d0
   tr_jbjb=0.0d0
   tr_jajb=0.0d0
   tr_jajc=0.0d0
   tr_jajd=0.0d0
   tr_jbjc=0.0d0
   tr_jbjd=0.0d0
!!!! what we want to compute is
   !!! <p1 p2| (Je_12)* | pp2 pp1>*< pp1 pp2| J_12 | p1 p2>
   do i1=1,2
      do f1=1,2
         J_1(f1,i1)=sum(ubarpp1(f1,:)*matmul(Pi_k1(:,:),up1(i1,:)))
         Je_1(f1,i1)=sum(ubarpp2(f1,:)*matmul(Pi_k1e(:,:),up1(i1,:)))
         Je_1_dag(f1,i1)=conjg(Je_1(f1,i1))
         J_2(f1,i1)=sum(ubarpp2(f1,:)*matmul(Pi_k2(:,:),up2(i1,:)))
         Je_2(f1,i1)=sum(ubarpp1(f1,:)*matmul(Pi_k2e(:,:),up2(i1,:)))
         Je_2_dag(f1,i1)=conjg(Je_2(f1,i1))
         do i=1,4
            J_12a(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_a_mu(:,:,i),up1(i1,:)))
            Je_12a(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(Je_a_mu(:,:,i),up1(i1,:)))
            J_12b(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_b_mu(:,:,i),up1(i1,:)))
            Je_12b(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(Je_b_mu(:,:,i),up1(i1,:)))
            Je_12a_dag(f1,i1,i)=conjg(Je_12a(f1,i1,i))
            Je_12b_dag(f1,i1,i)=conjg(Je_12b(f1,i1,i))
           !!
           J_12c(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(J_c_mu(:,:,i),up2(i1,:)))
           J_12d(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(J_d_mu(:,:,i),up2(i1,:)))
         enddo
      enddo
   enddo
   !....note that I am only evaluating < J^i> <J^i> because we are interested in the em responses and
   ! only the diagonal part contributes. In the electroweak case I will have 2 loops: i=1,4: j=1,4
   do i=1,4
      do j=1,4
         do i1=1,2
            do f1=1,2
               do i2=1,2
                  do f2=1,2
                     tr_jaja(i,j)=tr_jaja(i,j)+ Je_12a_dag(f2,i1,i)*J_12a(f1,i1,j)*Je_2_dag(f1,i2)*J_2(f2,i2)
                     tr_jbjb(i,j)=tr_jbjb(i,j)+ Je_12b_dag(f2,i1,i)*J_12b(f1,i1,j)*Je_2_dag(f1,i2)*J_2(f2,i2)
                     tr_jajb(i,j)=tr_jajb(i,j)+ Je_12a_dag(f2,i1,i)*J_12b(f1,i1,j)*Je_2_dag(f1,i2)*J_2(f2,i2)
                     tr_jajc(i,j)=tr_jajc(i,j)+ Je_12a_dag(f2,i1,i)*J_1(f1,i1)*Je_2_dag(f1,i2)*J_12c(f2,i2,j)
                     tr_jajd(i,j)=tr_jajd(i,j)+ Je_12a_dag(f2,i1,i)*J_1(f1,i1)*Je_2_dag(f1,i2)*J_12d(f2,i2,j)
                     tr_jbjc(i,j)=tr_jbjc(i,j)+ Je_12b_dag(f2,i1,i)*J_1(f1,i1)*Je_2_dag(f1,i2)*J_12c(f2,i2,j)
                     tr_jbjd(i,j)=tr_jbjd(i,j)+ Je_12b_dag(f2,i1,i)*J_1(f1,i1)*Je_2_dag(f1,i2)*J_12d(f2,i2,j)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   
   jaja(:,:)=tr_jaja(:,:)*16.0d0/3.0d0
   jbjb(:,:)=tr_jbjb(:,:)*(-16.0d0/9.0d0)
   jajb(:,:)=2.0d0*tr_jajb(:,:)*(16.0d0/9.0d0)
   jajc(:,:)=tr_jajc(:,:)*(-16.0d0/9.0d0)
   jajd(:,:)=tr_jajd(:,:)*(16.0d0/9.0d0)
   jbjc(:,:)=tr_jbjc(:,:)*(16.0d0/9.0d0)
   jbjd(:,:)=tr_jbjd(:,:)*(16.0d0/3.0d0)

   r_exc=jaja+jbjb+jajb+jajc+jajd+jbjc+jbjd
   r_cc=r_exc(1,1)
   r_cl=-0.5d0*(r_exc(1,4)+r_exc(4,1))
   r_ll=r_exc(4,4)
   r_t=r_exc(2,2)+r_exc(3,3)
   r_tp=-0.5d0*ci*(r_exc(2,3)-r_exc(3,2)) 
  
  return
end subroutine det_JaJc_exc





 subroutine det_JpiJaJb(r_cc,r_cl,r_ll,r_t,r_tp)
   implicit none
   integer*4 :: i1,f1,i2,f2,i,j
   real*8 :: r_cc,r_cl,r_ll,r_t,r_tp
   complex*16 :: J_12a(2,2,4),J_12b(2,2,4)   
   complex*16 :: J_12a_dag(2,2,4),J_12b_dag(2,2,4)
   complex*16 :: J_2(2,2),J_2_dag(2,2)
   complex*16 :: J_1(2,2),J_1_dag(2,2),J_s2(2,2,4),J_s2_dag(2,2,4)
   complex*16 :: J_f(2,2,4),J_s1(2,2,4),J_f_dag(2,2,4),J_s1_dag(2,2,4)
   complex*16 :: J_p2(2,2,4),J_p2_dag(2,2,4),J_p1(2,2,4),J_p1_dag(2,2,4)
   complex*16 :: jfja(4,4),jfjb(4,4),jsja(4,4),jsjb(4,4),jpja(4,4),jpjb(4,4),r_av(4,4)

   do i1=1,2
      do f1=1,2
         J_1(f1,i1)=sum(ubarpp1(f1,:)*matmul(Pi_k1(:,:),up1(i1,:)))
         J_1_dag(f1,i1)=conjg(J_1(f1,i1))
         J_2(f1,i1)=sum(ubarpp2(f1,:)*matmul(Pi_k2(:,:),up2(i1,:)))
         J_2_dag(f1,i1)=conjg(J_2(f1,i1))
         do i=1,4
           J_12a(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_a_mu(:,:,i),up1(i1,:)))
           J_12b(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_b_mu(:,:,i),up1(i1,:)))
           J_12a_dag(f1,i1,i)=conjg(J_12a(f1,i1,i))
           J_12b_dag(f1,i1,i)=conjg(J_12b(f1,i1,i))
           J_f(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_pif(:,:,i),up1(i1,:)))
           J_f_dag(f1,i1,i)=conjg(J_f(f1,i1,i))
           J_s1(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_sea1(:,:,i),up1(i1,:)))
           J_s1_dag(f1,i1,i)=conjg(J_s1(f1,i1,i))
           J_s2(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(J_sea2(:,:,i),up2(i1,:)))
           J_s2_dag(f1,i1,i)=conjg(J_s2(f1,i1,i))
           J_p1(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_pl1(:,:,i),up1(i1,:)))
           J_p1_dag(f1,i1,i)=conjg(J_p1(f1,i1,i))
           J_p2(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(J_pl2(:,:,i),up2(i1,:)))
           J_p2_dag(f1,i1,i)=conjg(J_p2(f1,i1,i))
         enddo
      enddo
   enddo   
   jfja=czero
   jfjb=czero
   jsja=czero
   jsjb=czero
   jpja=czero
   jpjb=czero
   do i=1,4
      do j=1,4
         do i1=1,2
            do f1=1,2
               do i2=1,2
                  do f2=1,2
                     jfja(i,j)=jfja(i,j)+J_12a_dag(f1,i1,i)*J_f(f1,i1,j)*J_2_dag(f2,i2)*J_2(f2,i2)
                     jfjb(i,j)=jfjb(i,j)+J_12b_dag(f1,i1,i)*J_f(f1,i1,j)*J_2_dag(f2,i2)*J_2(f2,i2)
                     jsja(i,j)=jsja(i,j)+J_12a_dag(f1,i1,i)*J_s1(f1,i1,j)*J_2_dag(f2,i2)*J_2(f2,i2) &
                          &  +J_12a_dag(f1,i1,i)*J_1(f1,i1)*J_2_dag(f2,i2)*J_s2(f2,i2,j)
                     jsjb(i,j)=jsjb(i,j)+J_12b_dag(f1,i1,i)*J_s1(f1,i1,j)*J_2_dag(f2,i2)*J_2(f2,i2) &
                          &  +J_12b_dag(f1,i1,i)*J_1(f1,i1)*J_2_dag(f2,i2)*J_s2(f2,i2,j)
                     jpja(i,j)=jpja(i,j)+J_12a_dag(f1,i1,i)*J_p1(f1,i1,j)*J_2_dag(f2,i2)*J_2(f2,i2) &
                          &  +J_12a_dag(f1,i1,i)*J_1(f1,i1)*J_2_dag(f2,i2)*J_p2(f2,i2,j)
                     jpjb(i,j)=jpjb(i,j)+J_12b_dag(f1,i1,i)*J_p1(f1,i1,j)*J_2_dag(f2,i2)*J_2(f2,i2) &
                          &  +J_12b_dag(f1,i1,i)*J_1(f1,i1)*J_2_dag(f2,i2)*J_p2(f2,i2,j)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   
   jfja=jfja*16.0d0/3.0d0*2.0d0
   jfjb=jfjb*(-16.0d0/3.0d0)*2.0d0
   jsja=jsja*16.0d0/3.0d0*2.0d0
   jsjb=jsjb*(-16.0d0/3.0d0)*2.0d0
   jpja=jpja*16.0d0/3.0d0*2.0d0
   jpjb=jpjb*(-16.0d0/3.0d0)*2.0d0

   r_av=jfja+jfjb+jsja+jsjb+jpja+jpjb
   r_cc=r_av(1,1)
   r_cl=-0.5d0*(r_av(1,4)+r_av(4,1))
   r_ll=r_av(4,4)
   r_t=r_av(2,2)+r_av(3,3)
   r_tp=-0.5d0*ci*(r_av(2,3)-r_av(3,2)) 

   return
 end subroutine det_JpiJaJb



   subroutine det_JpiJaJb_exc(r_cc,r_cl,r_ll,r_t,r_tp)
   implicit none
   integer*4 :: i1,f1,i2,f2,i,j
   real*8 :: r_cc,r_cl,r_ll,r_t,r_tp
   complex*16 :: J_12a(2,2,4),J_12b(2,2,4)
   complex*16 :: Je_12a(2,2,4),Je_12b(2,2,4)   
   complex*16 :: Je_12a_dag(2,2,4),Je_12b_dag(2,2,4)
   complex*16 :: J_2(2,2),Je_2(2,2),Je_2_dag(2,2)
   complex*16 :: J_1(2,2),J_1_dag(2,2),J_s2(2,2,4),J_s2_dag(2,2,4)
   complex*16 :: J_f(2,2,4),J_s1(2,2,4),J_f_dag(2,2,4),J_s1_dag(2,2,4)
   complex*16 :: J_p2(2,2,4),J_p2_dag(2,2,4),J_p1(2,2,4),J_p1_dag(2,2,4)
   complex*16 :: tr_2,tr_1
   complex*16 :: jfja(4,4),jfjb(4,4),jsja(4,4),jsjb(4,4),jpja(4,4),jpjb(4,4),r_av(4,4)   

   do i1=1,2
      do f1=1,2
         J_1(f1,i1)=sum(ubarpp1(f1,:)*matmul(Pi_k1(:,:),up1(i1,:)))
         J_1_dag(f1,i1)=conjg(J_1(f1,i1))
         J_2(f1,i1)=sum(ubarpp2(f1,:)*matmul(Pi_k2(:,:),up2(i1,:)))
         Je_2(f1,i1)=sum(ubarpp1(f1,:)*matmul(Pi_k2e(:,:),up2(i1,:)))
         Je_2_dag(f1,i1)=conjg(Je_2(f1,i1))
         do i=1,4
            J_12a(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_a_mu(:,:,i),up1(i1,:)))
            Je_12a(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(Je_a_mu(:,:,i),up1(i1,:)))
            Je_12a_dag(f1,i1,i)=conjg(Je_12a(f1,i1,i))
            J_12b(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_b_mu(:,:,i),up1(i1,:)))
            Je_12b(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(Je_b_mu(:,:,i),up1(i1,:)))
            Je_12b_dag(f1,i1,i)=conjg(Je_12b(f1,i1,i))
            J_f(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_pif(:,:,i),up1(i1,:)))
            J_f_dag(f1,i1,i)=conjg(J_f(f1,i1,i))
            J_s1(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_sea1(:,:,i),up1(i1,:)))
            J_s1_dag(f1,i1,i)=conjg(J_s1(f1,i1,i))
            J_s2(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(J_sea2(:,:,i),up2(i1,:)))
            J_s2_dag(f1,i1,i)=conjg(J_s2(f1,i1,i))
            J_p1(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_pl1(:,:,i),up1(i1,:)))
            J_p1_dag(f1,i1,i)=conjg(J_p1(f1,i1,i))
            J_p2(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(J_pl2(:,:,i),up2(i1,:)))
            J_p2_dag(f1,i1,i)=conjg(J_p2(f1,i1,i))
         enddo
      enddo
   enddo   
   jfja=czero
   jfjb=czero
   jsja=czero
   jsjb=czero
   jpja=czero
   jpjb=czero
   do i=1,4
      do j=1,4
         do i1=1,2
            do f1=1,2
               do i2=1,2
                  do f2=1,2
                     jfja(i,j)=jfja(i,j)+Je_12a_dag(f2,i1,i)*J_f(f1,i1,j)*Je_2_dag(f1,i2)*J_2(f2,i2)
                     jfjb(i,j)=jfjb(i,j)+Je_12b_dag(f2,i1,i)*J_f(f1,i1,j)*Je_2_dag(f1,i2)*J_2(f2,i2)
                     jsja(i,j)=jsja(i,j)+Je_12a_dag(f2,i1,i)*J_s1(f1,i1,j)*Je_2_dag(f1,i2)*J_2(f2,i2) &
                          &  +Je_12a_dag(f2,i1,i)*J_1(f1,i1)*Je_2_dag(f1,i2)*J_s2(f2,i2,j)
                     jsjb(i,j)=jsjb(i,j)+Je_12b_dag(f2,i1,i)*J_s1(f1,i1,j)*Je_2_dag(f1,i2)*J_2(f2,i2) &
                          &  +Je_12b_dag(f2,i1,i)*J_1(f1,i1)*Je_2_dag(f1,i2)*J_s2(f2,i2,j)
                     jpja(i,j)=jpja(i,j)+Je_12a_dag(f2,i1,i)*J_p1(f1,i1,j)*Je_2_dag(f1,i2)*J_2(f2,i2) &
                          &  +Je_12a_dag(f2,i1,i)*J_1(f1,i1)*Je_2_dag(f1,i2)*J_p2(f2,i2,j)
                     jpjb(i,j)=jpjb(i,j)+Je_12b_dag(f2,i1,i)*J_p1(f1,i1,j)*Je_2_dag(f1,i2)*J_2(f2,i2) &
                          &  +Je_12b_dag(f2,i1,i)*J_1(f1,i1)*Je_2_dag(f1,i2)*J_p2(f2,i2,j)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   jfja(:,:)=jfja(:,:)*16.0d0/3.0d0*2.0d0
   jfjb(:,:)=jfjb(:,:)*16.0d0/3.0d0*2.0d0
   jsja(:,:)=jsja(:,:)*16.0d0/3.0d0*2.0d0
   jsjb(:,:)=jsjb(:,:)*16.0d0/3.0d0*2.0d0
   jpja(:,:)=jpja(:,:)*16.0d0/3.0d0*2.0d0
   jpjb(:,:)=jpjb(:,:)*16.0d0/3.0d0*2.0d0
   
   r_av=jfja+jfjb+jsja+jsjb+jpja+jpjb
   r_cc=r_av(1,1)
   r_cl=-0.5d0*(r_av(1,4)+r_av(4,1))
   r_ll=r_av(4,4)
   r_t=r_av(2,2)+r_av(3,3)
   r_tp=-0.5d0*ci*(r_av(2,3)-r_av(3,2)) 

   return
 end subroutine det_JpiJaJb_exc

subroutine det_JpiJpi(r_cc,r_cl,r_ll,r_t,r_tp)
   implicit none
   integer*4 :: i1,f1,i2,f2,mu,nu
   real*8 :: r_cc,r_cl,r_ll,r_t,r_tp   
   complex*16 :: jfjf(4,4), jsjs(4,4),jfjs(4,4),jpjs(4,4),r_av(4,4) 
   complex*16 ::  jpjp(4,4),jfjp(4,4)
   complex*16 :: tr_2,tr_1   
   complex*16:: tr_ff(4,4),tr_ss11(4,4),tr_ss12(4)
   complex*16 :: tr_ss21(4),tr_ss22(4,4),tr_fs1(4,4)
   complex*16 :: tr_pp11(4,4),tr_pp12(4),tr_pp21(4),tr_pp22(4,4),tr_fp1(4,4)
   complex*16 :: tr_ps11(4,4),tr_ps12(4),tr_ps21(4),tr_ps22(4,4)
   complex*16 :: J_2(2,2),J_2_dag(2,2)
   complex*16 :: J_1(2,2),J_1_dag(2,2),J_s2(2,2,4),J_s2_dag(2,2,4)
   complex*16 :: J_f(2,2,4),J_s1(2,2,4),J_f_dag(2,2,4),J_s1_dag(2,2,4)
   complex*16 :: J_p2(2,2,4),J_p2_dag(2,2,4),J_p1(2,2,4),J_p1_dag(2,2,4)


   J_2=czero
   J_2_dag=czero
   tr_2=czero
   J_s2=czero
   J_s2_dag=czero
   J_p2=czero
   J_p2_dag=czero
   tr_ss21=czero
   tr_ss22=czero
   tr_pp21=czero
   tr_pp22=czero
   tr_ps22=czero
   tr_ps21=czero

   do i2=1,2
      do f2=1,2
         J_2(f2,i2)=sum(ubarpp2(f2,:)*matmul(Pi_k2(:,:),up2(i2,:)))
         J_2_dag(f2,i2)=conjg(J_2(f2,i2))
         tr_2=tr_2+J_2_dag(f2,i2)*J_2(f2,i2)
         do mu=1,4
           J_s2(f2,i2,mu)=sum(ubarpp2(f2,:)*matmul(J_sea2(:,:,mu),up2(i2,:)))
           J_s2_dag(f2,i2,mu)=conjg(J_s2(f2,i2,mu))
           J_p2(f2,i2,mu)=sum(ubarpp2(f2,:)*matmul(J_pl2(:,:,mu),up2(i2,:)))
           J_p2_dag(f2,i2,mu)=conjg(J_p2(f2,i2,mu))
         enddo
      enddo
   enddo   
   do i2=1,2
      do f2=1,2
         do mu=1,4
            tr_ss21(mu)=tr_ss21(mu)+J_s2_dag(f2,i2,mu)*J_2(f2,i2)
            tr_pp21(mu)=tr_pp21(mu)+J_p2_dag(f2,i2,mu)*J_2(f2,i2)
          !  tr_pp21(mu)=tr_pp21(mu)+J_p2_dag(f2,i2,mu)*J_2(f2,i2)
          !  tr_pp21(mu)=tr_pp21(mu)+J_2_dag(f2,i2)*J_s2(f2,i2,mu)
            do nu=1,4
               tr_ss22(mu,nu)=tr_ss22(mu,nu)+J_s2_dag(f2,i2,mu)*J_s2(f2,i2,nu)
               tr_pp22(mu,nu)=tr_pp22(mu,nu)+J_p2_dag(f2,i2,mu)*J_p2(f2,i2,nu)
               tr_ps22(mu,nu)=tr_ps22(mu,nu)+J_p2_dag(f2,i2,mu)*J_s2(f2,i2,nu)
            enddo
         enddo
      enddo 
   enddo  

   J_1=czero
   J_1_dag=czero
   tr_1=czero
   J_f=czero   
   J_f_dag=czero
   J_s1=czero
   J_s1_dag=czero
   J_p1=czero
   J_p1_dag=czero
   tr_ff=czero
   tr_ss11=czero
   tr_ss12=czero
   tr_fs1=czero
   tr_pp11=czero
   tr_pp12=czero
   tr_fp1=czero
   tr_ps12=czero
   tr_ps11=czero

   do i1=1,2
      do f1=1,2
         J_1(f1,i1)=sum(ubarpp1(f1,:)*matmul(Pi_k1(:,:),up1(i1,:)))
         J_1_dag(f1,i1)=conjg(J_1(f1,i1))
         tr_1=tr_1+J_1_dag(f1,i1)*J_1(f1,i1)
         do mu=1,4
         J_f(f1,i1,mu)=sum(ubarpp1(f1,:)*matmul(J_pif(:,:,mu),up1(i1,:)))
         J_s1(f1,i1,mu)=sum(ubarpp1(f1,:)*matmul(J_sea1(:,:,mu),up1(i1,:)))
         J_p1(f1,i1,mu)=sum(ubarpp1(f1,:)*matmul(J_pl1(:,:,mu),up1(i1,:)))
         J_f_dag(f1,i1,mu)=conjg(J_f(f1,i1,mu))
         J_s1_dag(f1,i1,mu)=conjg(J_s1(f1,i1,mu))
         J_p1_dag(f1,i1,mu)=conjg(J_p1(f1,i1,mu))
         enddo
      enddo  
   enddo 

   do i1=1,2
      do f1=1,2
         do mu=1,4
            tr_ss12(mu)=tr_ss12(mu)+J_1_dag(f1,i1)*J_s1(f1,i1,mu)
            tr_pp12(mu)=tr_pp12(mu)+J_1_dag(f1,i1)*J_p1(f1,i1,mu)
            tr_ps12(mu)=tr_ps12(mu)+J_1_dag(f1,i1)*J_s1(f1,i1,mu)
           ! tr_ps12(mu)=tr_ps12(mu)+J_s1_dag(f1,i1,mu)*J_1(f1,i1)
            do nu=1,4  
               tr_ff(mu,nu)=tr_ff(mu,nu)+J_f_dag(f1,i1,mu)*J_f(f1,i1,nu)
               tr_ss11(mu,nu)=tr_ss11(mu,nu)+J_s1_dag(f1,i1,mu)*J_s1(f1,i1,nu)
               tr_fs1(mu,nu)=tr_fs1(mu,nu)+J_f_dag(f1,i1,mu)*J_s1(f1,i1,nu)
               tr_pp11(mu,nu)=tr_pp11(mu,nu)+J_p1_dag(f1,i1,mu)*J_p1(f1,i1,nu)
               tr_ps11(mu,nu)=tr_ps11(mu,nu)+J_p1_dag(f1,i1,mu)*J_s1(f1,i1,nu)
               tr_fp1(mu,nu)=tr_fp1(mu,nu)+J_f_dag(f1,i1,mu)*J_p1(f1,i1,nu)
            enddo
         enddo   
      enddo
   enddo

   jfjf(:,:)=16.0d0*tr_ff(:,:)*tr_2
   do mu=1,4
      do nu=1,4
      jsjs(mu,nu)=2.0d0*tr_ss12(mu)*tr_ss21(nu)
      jpjp(mu,nu)=2.0d0*tr_pp12(mu)*tr_pp21(nu)
      jpjs(mu,nu)=2.0d0*tr_ps12(mu)*tr_ps21(nu)
      enddo
   enddo 
   jsjs(:,:)=16.0d0*(jsjs(:,:)+tr_ss11(:,:)*tr_2+tr_ss22(:,:)*tr_1)
   jpjp(:,:)=16.0d0*(jpjp(:,:)+tr_pp11(:,:)*tr_2+tr_pp22(:,:)*tr_1)
   jpjs(:,:)=16.0d0*(jpjs(:,:)+tr_ps11(:,:)*tr_2+tr_ps22(:,:)*tr_1)
   jfjs(:,:)=16.0d0*(tr_fs1(:,:)*tr_2*2.0d0)
   jfjp(:,:)=16.0d0*(tr_fp1(:,:)*tr_2*2.0d0)
   r_av=jfjf+jsjs+2.0d0*jfjs+jpjp+2.0d0*jfjp+2.0d0*jpjs
   
   r_cc=r_av(1,1)
   r_cl=-0.5d0*(r_av(1,4)+r_av(4,1))
   r_ll=r_av(4,4)
   r_t=r_av(2,2)+r_av(3,3)
   r_tp=-0.5d0*ci*(r_av(2,3)-r_av(3,2)) 

   return
end subroutine   


subroutine delta_se(pd2,width,pot)
   implicit none
   real*8 :: pd2,width,kpi,ekpi,eknuc,r2a,rfa,pot
   width=0.0d0
   !pot=0.0d0!-40.0d0
   if (pd2.ge.(xmpi+xmn)**2)then
      kpi=dsqrt(1.0d0/4.0d0/pd2*(pd2-(xmn+xmpi)**2)*(pd2-(xmn-xmpi)**2))
      ekpi=dsqrt(kpi**2+xmpi**2)
      eknuc=dsqrt(kpi**2+xmn**2)
      r2a= (eknuc-ekpi)**2 -4.0d0*kpi**2
      rfa=sqrt(0.95d0*xmn**2/(0.95d0*xmn**2-r2a))
!.....definizione Gamma
!      width=(4.0d0*fpind)**2/12.0d0/pi/xmpi**2*kpi**3/sqrt(pd2)*(xmn+eknuc)  &
!   &  *rfa**2!*(lpind**2/(lpind**2-xmpi**2))**2
!      width=120.0d0

      width=0.38/(3.0d0*xmpi**2)*kpi**3/sqrt(pd2)*(xmn+eknuc)  &
   &  *rfa**2-pot*2.0d0!*(lpind**2/(lpind**2-xmpi**2))**2


      return
   endif
!.....Eq.3.26 Phys.Rev. C 49 2650

   return
end subroutine


end module

