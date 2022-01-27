     subroutine nform(ichoice,Q2,F_1s,F_2s,F_1v,F_2v,G_es,G_ms,G_ev,G_mv)

!      
!     revised January 2014, Virginia Tech
!
!     ichoice = 1  Dipole parametrization
!     ichoice = 2  Kelly's parametrization [PRC 70, 068282 (2004)]
!     ichoice = 3  BBBA parametrization [NPB Proc. Suppl., 159, 127 (2006)]
!     ichoice =4   Arrington-Hill parametrization [Phys.Lett.B, 777, pgg. 8-15 (2018)]       
       
!
!     Q2 is given in units of fm^-2
!
     implicit none
     integer*4 :: ichoice,i,j
     real*8, parameter :: tcut = 0.0779191396,t0 = -0.7d0
     real*8 :: Q2,Q2G,G_ep,G_mp,G_en,G_mn,F_1s,F_1v,F_2s,F_2v,G_D &
    &       ,Lambdasq,hc,pm,zero,mu_p,mu_n,a_ep,a_mp,a_mn,tau &
    &       ,a_en,b_en,num_ep,den_ep,num_mp,den_mp,num_en,den_en &
    &       ,num_mn,den_mn,one,G_es,G_ev,G_ms,G_mv,G_A,gan1,MAsq
     real*8 :: b_ep(3),b_mp(3),b_mn(3),aa_ep(4),aa_mp(4),aa_en(4) &
          &      ,aa_mn(4),bb_ep(4),bb_mp(4),bb_en(4),bb_mn(4)
     real*8 :: gep_coef(13),gmp_coef(13),gen_coef(13),gmn_coef(13)
     real*8 :: z
!     
     data a_ep , a_mp ,a_mn  / -0.24 , 0.12 , 2.33 / 
     data b_ep / 10.98 , 12.82 , 21.97 /
     data b_mp / 10.97 , 18.86 , 6.55 /
     data b_mn / 14.72 , 24.20 , 84.1 /
     data a_en , b_en / 1.70 , 3.30 /
!
     data aa_ep / 1.0 , -0.0578 , 0.0  , 0.0 /     
     data aa_mp / 1.0 ,  0.0150 , 0.0  , 0.0 /
     data aa_en / 0.0 ,  1.2500 , 1.30 , 0.0 /
     data aa_mn / 1.0 ,  1.8100 , 0.0  , 0.0 /
!
     data bb_ep / 11.10 ,  13.60 ,   33.00 ,   0.0 /
     data bb_mp / 11.10 ,  19.60 ,    7.54 ,   0.0 /
     data bb_en / 9.86  , 305.00 , -758.00 , 802.0 /
     data bb_mn / 14.10 ,  20.70 ,   68.70 ,   0.0 /
!     
     data gep_coef / 0.239163298067, -1.10985857441, 1.44438081306, 0.479569465603, &
          & -2.28689474187,  1.12663298498, 1.25061984354,-3.63102047159, 4.08221702379, &
          & 0.504097346499,  -5.08512046051,  3.96774254395,-0.981529071103 /
     data gmp_coef / 0.264142994136, -1.09530612212, 1.21855378178, 0.661136493537, -1.40567892503, &
          &  -1.35641843888, 1.44702915534, 4.2356697359, -5.33404565341, -2.91630052096,    &
          & 8.70740306757, -5.70699994375, 1.28081437589 /
     data gen_coef / 0.048919981379,-0.064525053912,-0.240825897382,0.392108744873, 0.300445258602, &
          & -0.661888687179,-0.175639769687, 0.624691724461,-0.077684299367,-0.236003975259, 0.090401973470, &
          0.0, 0.0 /
     data gmn_coef / 0.257758326959,-1.079540642058, 1.182183812195,0.711015085833,-1.348080936796,-1.662444025208, &
          & 2.624354426029, 1.751234494568,-4.922300878888, 3.197892727312,-0.712072389946, 0.0, 0.0 /
!     
     hc = .197327d0
     pm = .938272d0    ! proton mass in GeV
     zero = 0.d0
     one = 1.d0
     mu_p = 2.79278d0   ! proton magnetic moment
     mu_n = -1.91315d0  ! neutron magnetic moment
     Q2G = Q2*hc*hc
     tau = Q2G/4./pm/pm
     gan1= 1.2694d0
!
!    Omar's value 
!    Lambdasq = 0.71
!    Rocco's value (I believe from Ciofi)
!     Lambdasq = 0.69513d0
     Lambdasq = 0.7174d0
     MAsq=1.0d0

!     Lambdasq = 0.55d0

     G_D=1.0d0/(1.0d0+Q2G/Lambdasq)**2
     G_A=gan1/(1.0d0+Q2G/MAsq)**2
!    
     if(ichoice.eq.1)then
!             
       G_ep = G_D 
!       G_en = zero
       G_en = -mu_n*Q2G*G_D/(1.+Q2G/pm**2)/(4.*pm**2)

       G_mp = mu_p*G_D
       G_mn = mu_n*G_D
     else if(ichoice.eq.2)then
       G_ep = (1. + a_ep*tau)/ &
    &         (1. + b_ep(1)*tau + b_ep(2)*tau**2 + b_ep(3)*tau**3) 
       G_en = G_D*a_en*tau/(1. + b_en*tau)
       G_mp = mu_p*(1. + a_mp*tau)/ &
    &              (1. + b_mp(1)*tau + b_mp(2)*tau**2 + b_mp(3)*tau**3)
       G_mn = mu_n*(1. + a_mn*tau)/ &
    &              (1. + b_mn(1)*tau + b_mn(2)*tau**2 + b_mn(3)*tau**3)
     else if(ichoice.eq.3)then   
       num_ep = zero
       den_ep = one
       num_mp = zero 
       den_mp = one
       num_en = zero 
       den_en = one
       num_mn = zero
       den_mn = one
       do i = 1,4
       j = i-1
       num_ep = num_ep + aa_ep(i)*tau**j
       den_ep = den_ep + bb_ep(i)*tau**i
       num_mp = num_mp + aa_mp(i)*tau**j
       den_mp = den_mp + bb_mp(i)*tau**i
       num_en = num_en + aa_en(i)*tau**j
       den_en = den_en + bb_en(i)*tau**i
       num_mn = num_mn + aa_mn(i)*tau**j
       den_mn = den_mn + bb_mn(i)*tau**i
       end do
       G_ep = num_ep/den_ep
       G_mp = mu_p*num_mp/den_mp
       G_en = num_en/den_en
       G_mn = mu_n*num_mn/den_mn
    else if(ichoice.eq.4) then
       G_ep=0.0d0
       G_en=0.0d0
       G_mp=0.0d0
       G_mn=0.0d0
       z = (sqrt(tcut+Q2G)-sqrt(tcut-t0))/(sqrt(tcut+Q2G)+sqrt(tcut-t0)) 
       do i=1,13
          G_ep=G_ep+gep_coef(i)*z**(i-1)
          G_mp=G_mp+gmp_coef(i)*z**(i-1)*mu_p
          G_en=G_en+ gen_coef(i)*z**(i-1)
          G_mn=G_mn+gmn_coef(i)*z**(i-1)*mu_n
       enddo
     end if
       G_es=G_ep+G_en
       G_ev=G_ep-G_en
       G_ms=G_mp+G_mn
       G_mv=G_mp-G_mn
       F_1s = (G_es+tau*G_ms)/(1.+tau)
       F_1v = (G_ev+tau*G_mv)/(1.+tau)
       F_2s = (G_ms-G_es)/(1.+tau)
       F_2v = (G_mv-G_ev)/(1.+tau)

!
     return
     end   
