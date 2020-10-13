module dirac_matrices
    implicit none
    integer*4, private, save :: i_fl
    complex*16, private, parameter :: czero = (0.0d0,0.0d0)
    complex*16, private, parameter :: cone  = (1.0d0,0.0d0)
    complex*16, private, parameter :: ci    = (0.0d0,1.0d0)
    real*8, private, parameter :: pi=acos(-1.0d0)    
    real*8, private, save :: mqe, qval
    complex*16, private, save :: sig(3,2,2),id(2,2),id4(4,4),up(2),down(2)
    complex*16, private, save :: up1(2,4),upp1(2,4), &
            &   ubarp1(2,4),ubarpp1(2,4)
    complex*16, private, save :: gamma_mu(4,4,5),g_munu(4,4),sigma_munu(4,4,4,4)
    complex*16, private, save :: q_sl(4,4)
    real*8, private, save ::  p1(4),pp1(4),q(4)
    complex*16, private, save :: J_1(4,4,4)
    real*8, private,save :: xmn
contains

subroutine dirac_matrices_in(xmn_in)
    implicit none
    integer*4 :: i,j
    real*8 :: xmd_in,xmn_in,xmpi_in
    xmn=xmn_in
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
    do i=1,4
       do j=1,4
          sigma_munu(:,:,i,j)=ci*0.5d0*(matmul(gamma_mu(:,:,i),gamma_mu(:,:,j)) &
               &     -matmul(gamma_mu(:,:,j),gamma_mu(:,:,i)))
       enddo
    enddo
          
    up(1)=cone;up(2)=czero
    down(1)=czero;down(2)=cone
end subroutine 

subroutine define_spinors()
    implicit none
    integer*4 :: i
    complex*16 :: sigp1(2,2),sigp2(2,2),sigpp1(2,2),sigpp2(2,2)
    real*8 :: cp1,cp2,cpp1,cpp2
    sigp1=czero
    sigpp1=czero
    !.....initialize quadrispinors
    up1=czero
    upp1=czero
!.......initialize normalization factors
    cp1=sqrt((p1(1)+xmn)/(2.0d0*p1(1)))
    cpp1=sqrt((pp1(1)+xmn)/(2.0d0*pp1(1)))
!.....define sigma*p
    do i=1,3
      sigp1=sigp1+sig(i,:,:)*p1(i+1)
      sigpp1=sigpp1+sig(i,:,:)*pp1(i+1)
    enddo
!.....build quadri-spinors    
    up1(1,1:2)=up(:)
    up1(1,3:4)=matmul(sigp1(:,:),up(:))/(p1(1)+xmn)
    up1(2,1:2)=down(:)
    up1(2,3:4)=matmul(sigp1(:,:),down(:))/(p1(1)+xmn)
    up1(:,:)=cp1*up1(:,:)
!
    upp1(1,1:2)=up(:)
    upp1(1,3:4)=matmul(sigpp1(:,:),up(:))/(pp1(1)+xmn)
    upp1(2,1:2)=down(:)
    upp1(2,3:4)=matmul(sigpp1(:,:),down(:))/(pp1(1)+xmn)
    upp1(:,:)=cpp1*upp1(:,:)
!
    ubarp1(1,1:2)=up(:)
    ubarp1(1,3:4)=-matmul(up(:),sigp1(:,:))/(p1(1)+xmn)
    ubarp1(2,1:2)=down(:)
    ubarp1(2,3:4)=-matmul(down(:),sigp1(:,:))/(p1(1)+xmn)
    ubarp1(:,:)=cp1*ubarp1(:,:)
!
    ubarpp1(1,1:2)=up(:)
    ubarpp1(1,3:4)=-matmul(up(:),sigpp1(:,:))/(pp1(1)+xmn)
    ubarpp1(2,1:2)=down(:)
    ubarpp1(2,3:4)=-matmul(down(:),sigpp1(:,:))/(pp1(1)+xmn)
    ubarpp1(:,:)=cpp1*ubarpp1(:,:)
    return
end subroutine

subroutine current_init(p1_in,pp1_in,q_in)
    implicit none
    real*8 ::  p1_in(4),pp1_in(4),q_in(4)
    p1=p1_in
    pp1=pp1_in
    q=q_in
    return
end subroutine




subroutine det_Ja(f1v,f2v)
  implicit none
  integer*4 :: mu,nu
  real*8 :: f1v,f2v  
  
  do mu=1,4
     J_1(:,:,mu)=czero
     do nu=1,4
        J_1(:,:,mu)=J_1(:,:,mu)+ci*f2v*sigma_munu(:,:,mu,nu)&
             &     *g_munu(nu,nu)*q(nu)/2.0d0/xmn
     enddo
     J_1(:,:,mu)=J_1(:,:,mu)+f1v*gamma_mu(:,:,mu)
  enddo
  
end subroutine det_Ja


subroutine det_res1b(rl,rt)
   implicit none
   integer*4 :: i1,f1,i,j
   complex*16 :: J_mu(2,2,4),J_mu_dag(2,2,4)
   real*8 :: res(4,4),rt,rl
   

   do i1=1,2
      do f1=1,2
         do i=1,4
            J_mu(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_1(:,:,i),up1(i1,:)))
            J_mu_dag(f1,i1,i)=conjg(J_mu(f1,i1,i))
         enddo
      enddo
   enddo
   
   res=0.0d0
   do i1=1,2
      do f1=1,2
         do i=1,4
            res(i,i)=res(i,i)+J_mu_dag(f1,i1,i)*J_mu(f1,i1,i)
         enddo
      enddo
   enddo
   
   rl=res(1,1)
   rt=res(2,2)+res(3,3)
 
  return
end subroutine det_res1b
end module dirac_matrices





