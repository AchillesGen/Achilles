module dirac_matrices_pi
    implicit none
    integer*4, private, parameter :: nspin_in=2,nspin_f=2    
    integer*4, private,parameter :: id_p=1,id_n=-1,id_pip=1,id_pi0=0,id_pim=-1   
    real*8, private,save :: xmn,xmp,xm,xmeta,xmpi0,xmpi1,pi,hbarc,xmnuc
    real*8, private, save ::  p1(4),pp1(4),q(4),kpi(4)
contains

subroutine dirac_matrices_in(xmp_in,xmn_in,xmeta_in,xmpi0_in,xmpi1_in,pi_in,hbarc_in)
    implicit none
    real*8 :: xmp_in,xmn_in,xmeta_in,xmpi0_in,xmpi1_in,pi_in,hbarc_in

    xmp=xmp_in
    xmn=xmn_in
    xmeta=xmeta_in
    xmpi0=xmpi0_in
    xmpi1=xmpi1_in
    pi=pi_in
    hbarc=hbarc_in
    xmnuc = (xmn + xmp)/2.0d0
end subroutine 

subroutine current_init(p1_in,pp1_in,q_in,kpi_in)
    implicit none
    real*8 ::  p1_in(4),pp1_in(4),q_in(4),kpi_in(4)
    real*8 :: w
    p1=p1_in
    pp1=pp1_in
    q=q_in
    w=q(1)
    kpi=kpi_in
    q(1)=w+p1(1)
    p1(1)=sqrt(p1(2)**2+p1(3)**2+p1(4)**2+xmnuc**2) 
    q(1)=q(1)-p1(1)
    
    return
end subroutine

subroutine hadr_curr_matrix_el(hid1,hid2,mesid1,has_axial,J_mu)
   implicit none     
   integer*4 :: i1,f1,i,j,ip
   integer*8 :: hid1,hid2,mesid1,DCC_mode,sumIDS
   real*8 :: PpiN(4), piNinv, Q2
   complex*16 :: zjmup(-1:1,-1:1,4)
   complex*16 :: J_mu(nspin_f,nspin_in,4)
   logical :: has_axial

   ! print*, mesid1, hid1, hid2
   if(hid1.eq.2212) then
           hid1=id_p
   elseif(hid1.eq.2112) then
           hid1=id_n
   endif

   if(hid2.eq.2212) then
           hid2=id_p
   elseif(hid2.eq.2112) then
           hid2=id_n
   endif

   if(mesid1.eq.111) then
           mesid1=id_pi0
   elseif(mesid1.eq.211) then
           mesid1=id_pip
   elseif(mesid1.eq.-211) then
           mesid1=id_pim       
   endif

   ! Sum the hadron id's to decide 
   ! which process we are doing
   sumIDS = hid1 + hid2 + mesid1
   DCC_mode = 0

   ! Decide which mode we are using
   ! Is the res axial form factor = 0?
   if (has_axial .eqv. .false.) then
        DCC_mode = 10
   else if (abs(sumIDS).eq.3) then
        if (sumIDS.gt.0) then
                DCC_mode = 1
        else 
                DCC_mode = 3
        endif
   else if (sumIDS.eq.0) then
        if (hid1.lt.0) then
                DCC_mode = 1
        else 
                DCC_mode = 3
        endif
   else if (abs(sumIDS).eq.2) then
        DCC_mode = -1 
   else if (abs(sumIDS).eq.1.and.hid1.eq.hid2) then
        if (hid1.lt.0) then
                DCC_mode = 1 
        else
                DCC_mode = 3
        endif
   else
        DCC_mode = -1
   endif

   PpiN = kpi + pp1
   piNinv = PpiN(1)**2 - PpiN(2)**2 - PpiN(3)**2 - PpiN(4)**2
   piNinv = abs(piNinv)

   Q2 = q(2)**2 + q(3)**2 + q(4)**2 - q(1)**2
   if (Q2.lt.0.0d0) then        
        J_mu(:,:,:) = (0.0d0,0.0d0)
        return
   endif

   if(sqrt(piNinv).lt.1076.957d0.or.sqrt(piNinv).gt.1400.0d0) then
        J_mu(:,:,:) = (0.0d0,0.0d0)
        return
   endif
   ! print*, "DCC_mode = ", DCC_mode

   call amplitude(q,pp1,kpi,DCC_mode,0,hid1,mesid1,1,zjmup(:,:,:))

   J_mu(1,1,:)=zjmup(-1,-1,:)
   J_mu(1,2,:)=zjmup(-1,1,:)
   J_mu(2,1,:)=zjmup(1,-1,:)  
   J_mu(2,2,:)=zjmup(1,1,:)
   J_mu=J_mu/hbarc*xmnuc/(2.0d0*pi)**3
   return
end subroutine
end module dirac_matrices_pi
