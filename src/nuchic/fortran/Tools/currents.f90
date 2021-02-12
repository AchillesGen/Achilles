module dirac_matrices
    implicit none

    private
    public :: dirac_matrix
    complex*16, parameter :: czero = (0.0d0,0.0d0)
    complex*16, parameter :: cone  = (1.0d0,0.0d0)
    complex*16, parameter :: ci    = (0.0d0,1.0d0)
    real*8, parameter :: pi=acos(-1.0d0)    

    type dirac_matrix_factory
        private
            character(len=20) :: factory_type
            class(dirac_matrix), pointer :: dirac_matrix_type
        contains
            procedure :: init
            procedure :: create_dirac_matrix
            procedure :: final
    end type dirac_matrix_factory

    type, abstract :: dirac_matrix
        contains
            procedure(gen_initialize), deferred :: initialize
            procedure(gen_spinors), deferred :: spinors
            procedure(gen_current), deferred :: current
            procedure(gen_Ja), deferred :: Ja
            procedure(gen_res1b), deferred :: res1b
    end type dirac_matrix

    abstract interface
        subroutine gen_initialize(self, xmn_in, hbarc_in)
            import :: dirac_matrix
            implicit none
            class(dirac_matrix), intent(in) :: self
            real*8, intent(in) :: xmn_in, hbarc_in
        end subroutine

        subroutine gen_spinors(self)
            import :: dirac_matrix
            implicit none
            class(dirac_matrix), intent(in) :: self
        end subroutine

        subroutine gen_current(self, mom)
            import :: dirac_matrix
            implicit none
            class(dirac_matrix), intent(in) :: self
            real*8, intent(in) :: mom(:, :)
        end subroutine

        subroutine gen_ja(self, fv)
            import :: dirac_matrix
            implicit none
            class(dirac_matrix), intent(in) :: self
            real*8, intent(in) :: fv(*)
        end subroutine

        subroutine gen_res1b(self, rl, rt)
            import :: dirac_matrix
            implicit none
            class(dirac_matrix), intent(in) :: self
            real*8, intent(in) :: rl(:,:), rt(:,:)
        end subroutine
    end interface

    type, extends(dirac_matrix) :: one_body
            complex*16 :: sig(3, 2, 2), id(2, 2), id4(4, 4), up(2), down(2)
            complex*16 :: up1(2, 4), upp1(2, 4), &
                & ubarp1(2, 4), ubarpp1(2, 4)
            complex*16 :: gamma_mu(4, 4, 5), g_munu(4, 4), sigma_munu(4, 4, 4, 4)
            complex*16 :: q_sl(4, 4), J_1(4, 4, 4)
            real*8 :: p1(4), pp1(4), q(4), xmn
        contains
            procedure :: initialize => onebody_init
            procedure :: spinors => onebody_spinors
            procedure :: current => onebody_current
            procedure :: Ja => onebody_ja
            procedure :: res1b => onebody_res1b
    end type one_body

    type, extends(dirac_matrix) :: pion_production
            real*8 :: xmn, hbarc
            real*8 :: p1(4),pp1(4),q(4),kpi(4)
        contains
            procedure :: initialize => pion_init
            procedure :: spinors => pion_spinors
            procedure :: current => pion_current
            procedure :: Ja => pion_ja
            procedure :: res1b => pion_res1b
    end type pion_production

contains

    subroutine init(self, string)
        class(dirac_matrix_factory), intent(inout) :: self
        character(len=*), intent(in) :: string
        self%factory_type = trim(string)
        self%dirac_matrix_type => null()
    end subroutine init

    subroutine final(self)
        class(dirac_matrix_factory), intent(inout) :: self
        deallocate(self%dirac_matrix_type)
        nullify(self%dirac_matrix_type)
    end subroutine final

    function create_dirac_matrix(self) result(ptr)
        class(dirac_matrix_factory) :: self
        class(dirac_matrix), pointer :: ptr

        if(self%factory_type == "OneBody") then
            if(associated(self%dirac_matrix_type)) deallocate(self%dirac_matrix_type)
            allocate(one_body :: self%dirac_matrix_type)
            ptr => self%dirac_matrix_type
        else if(self%factory_type == "PionProduction") then
            if(associated(self%dirac_matrix_type)) deallocate(self%dirac_matrix_type)
            allocate(pion_production :: self%dirac_matrix_type)
            ptr => self%dirac_matrix_type
        end if
    end function create_dirac_matrix

    subroutine onebody_init(self, xmn_in, hbarc_in)
        implicit none
        integer*4 :: i,j
        class(one_body), intent(inout) :: self
        real*8, intent(in) :: xmn_in,hbarc_in
        self%xmn=xmn_in
        self%sig(:,:,:)=czero
        self%id(:,:)=czero
        self%id(1,1)=cone;self%id(2,2)=cone
        self%sig(1,1,2)=cone;self%sig(1,2,1)=cone
        self%sig(2,1,2)=-ci;self%sig(2,2,1)=ci
        self%sig(3,1,1)=cone;self%sig(3,2,2)=-cone
        self%gamma_mu=czero
        self%gamma_mu(1:2,1:2,1)=self%id;self%gamma_mu(3:4,3:4,1)=-id
        self%id4=czero    
        self%id4(1:2,1:2)=self%id;self%id4(3:4,3:4)=self%id
        do i=2,4
          self%gamma_mu(1:2,3:4,i)=self%sig(i-1,:,:)
          self%gamma_mu(3:4,1:2,i)=-self%sig(i-1,:,:)
        enddo
        self%gamma_mu(1:2,3:4,5)=self%id
        self%gamma_mu(3:4,1:2,5)=self%id
        self%g_munu=czero
        self%g_munu(1,1)=cone;self%g_munu(2,2)=-cone;self%g_munu(3,3)=-cone;self%g_munu(4,4)=-cone
        do i=1,4
           do j=1,4
              self%sigma_munu(:,:,i,j)=ci*0.5d0*(matmul(self%gamma_mu(:,:,i),self%gamma_mu(:,:,j)) &
                   &     -matmul(self%gamma_mu(:,:,j),self%gamma_mu(:,:,i)))
           enddo
        enddo
              
        self%up(1)=cone;self%up(2)=czero
        self%down(1)=czero;self%down(2)=cone
    end subroutine 

    subroutine onebody_spinors(self)
        implicit none
        class(one_body), intent(inout) :: self
        integer*4 :: i
        complex*16 :: sigp1(2,2),sigp2(2,2),sigpp1(2,2),sigpp2(2,2)
        real*8 :: cp1,cp2,cpp1,cpp2
        sigp1=czero
        sigpp1=czero
        !.....initialize quadrispinors
        self%up1=czero
        self%upp1=czero
    !.......initialize normalization factors
        cp1=sqrt((self%p1(1)+self%xmn)/(2.0d0*self%p1(1)))
        cpp1=sqrt((self%pp1(1)+self%xmn)/(2.0d0*self%pp1(1)))
    !.....define sigma*p
        do i=1,3
          sigp1=sigp1+self%sig(i,:,:)*self%p1(i+1)
          sigpp1=sigpp1+self%sig(i,:,:)*self%pp1(i+1)
        enddo
    !.....build quadri-spinors    
        self%up1(1,1:2)=self%up(:)
        self%up1(1,3:4)=matmul(sigp1(:,:),self%up(:))/(self%p1(1)+self%xmn)
        self%up1(2,1:2)=self%down(:)
        self%up1(2,3:4)=matmul(sigp1(:,:),self%down(:))/(self%p1(1)+self%xmn)
        self%up1(:,:)=cp1*self%up1(:,:)
    !
        self%upp1(1,1:2)=self%up(:)
        self%upp1(1,3:4)=matmul(sigpp1(:,:),self%up(:))/(self%pp1(1)+self%xmn)
        self%upp1(2,1:2)=self%down(:)
        self%upp1(2,3:4)=matmul(sigpp1(:,:),self%down(:))/(self%pp1(1)+self%xmn)
        self%upp1(:,:)=cpp1*self%upp1(:,:)
    !
        self%ubarp1(1,1:2)=self%up(:)
        self%ubarp1(1,3:4)=-matmul(self%up(:),sigp1(:,:))/(self%p1(1)+self%xmn)
        self%ubarp1(2,1:2)=self%down(:)
        self%ubarp1(2,3:4)=-matmul(self%down(:),sigp1(:,:))/(self%p1(1)+self%xmn)
        self%ubarp1(:,:)=cp1*self%ubarp1(:,:)
    !
        self%ubarpp1(1,1:2)=self%up(:)
        self%ubarpp1(1,3:4)=-matmul(self%up(:),sigpp1(:,:))/(self%pp1(1)+self%xmn)
        self%ubarpp1(2,1:2)=self%down(:)
        self%ubarpp1(2,3:4)=-matmul(self%down(:),sigpp1(:,:))/(self%pp1(1)+self%xmn)
        self%ubarpp1(:,:)=cpp1*self%ubarpp1(:,:)
    end subroutine

    subroutine onebody_current(self, mom)
        implicit none
        class(one_body), intent(in) :: self
        real*8, intent(in) :: mom(:, :)
        p1 = mom(1)
        pp1 = mom(2)
        q = mom(3)
    end subroutine

    subroutine onebody_ja(self, fv)
        implicit none
        class(one_body), intent(in) :: self
        real*8, intent(in) :: fv(:)
        integer*4 :: mu, nu
        do mu=1,4
            J_1(:,:,mu)=czero
            do nu=1,4
                J_1(:,:,mu)=J_1(:,:,mu)+ci*fv(2)*sigma_munu(:,:,mu,nu)&
                    &     *g_munu(nu,nu)*q(nu)/2.0d0/xmn
            enddo
            J_1(:,:,mu)=J_1(:,:,mu)+fv(1)*gamma_mu(:,:,mu)
        enddo
    end subroutine onebody_ja

    subroutine onebody_res1b(self, rl, rt)
        implicit none
        class(one_body), intent(in) :: self
        real*8, intent(in) :: rl(2, 1), rt(2, 1)
        integer*4 :: i1,f1,i,j
        complex*16 :: J_mu(2,2,4),J_mu_dag(2,2,4)
        real*8 :: res(4,4)

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
        
        rl(1, 1)=rl(2, 1)=res(1,1)
        rt(1, 1)=rt(2, 1)=res(2,2)+res(3,3)
    end subroutine onebody_res1b

    subroutine pion_init(self, xmn_in, hbarc_in)
        implicit none
        class(pion_production), intent(in) :: self
        real*8, intent(in) :: xmn_in, hbarc_in
        xmn=xmn_in
        hbarc=hbarc_in
    end subroutine 

    subroutine pion_spinors(self)
        implicit none
        class(pion_production), intent(in) :: self
    end subroutine

    subroutine pion_current(self, mom)
        implicit none
        class(pion_production), intent(in) :: self
        real*8, intent(in) :: mom(4, 4)
        p1=mom(1)
        pp1=mom(2)
        q=mom(3)
        kpi=mom(4)
    end subroutine

    subroutine pion_ja(self, fv)
        implicit none
        class(pion_production), intent(in) :: self
        real*8, intent(in) :: fv(2)
    end subroutine

    subroutine pion_res1b(self, rl, rt)
        implicit none
        class(pion_production), intent(in) :: self
        real*8, intent(in) :: rt(2, 6),rl(2, 6)
        integer*4 :: i1,f1,i,j,ip
        complex*16 :: J_mu_p(3,6,3,4),J_mu_p_dag(3,6,3,4)
        complex*16 :: J_mu_n(3,6,3,4),J_mu_n_dag(3,6,3,4)
        real*8 :: res_p(4,4,6),res_n(4,4,6)
    ! n + pi- = 1
    ! n + pi0 = 2
    ! n + pi+ = 3
    ! p + pi- = 4
    ! p + pi0 = 5
    ! p + pi+ = 6
    
    !  call mxdelta(kpi,pp1,p1,q,J_mu_p,J_mu_n)
       !   call  mxgnpin(kpi,pp1,p1,q,J_mu_p,J_mu_n)
       !   call slmodel(kpi,pp1,p1,q,J_mu_p,J_mu_n)
       J_mu_p=(0.0d0,0.0d0)
       J_mu_n=(0.0d0,0.0d0)
       call amplitude(q,pp1,kpi,10,0,1,0,1,J_mu_p(:,3,:,:))
       call amplitude(q,pp1,kpi,10,0,1,1,1,J_mu_p(:,6,:,:))
       call amplitude(q,pp1,kpi,10,0,-1,0,1,J_mu_n(:,2,:,:))
       call amplitude(q,pp1,kpi,10,0,-1,-1,1,J_mu_n(:,4,:,:))
       J_mu_p=J_mu_p/hbarc
       J_mu_n=J_mu_n/hbarc
       
     J_mu_p_dag=0.0d0  
     do i1=1,3
        do ip=1,6
           do f1=1,3
              do i=1,4
                 J_mu_p_dag(f1,ip,i1,i)=conjg(J_mu_p(f1,ip,i1,i))
                ! if(ip.eq.1) write(100,*) 'jmup',J_mu_p_dag(f1,1,i1,i),ip
                ! if(ip.eq.2) write(101,*) 'jmup',J_mu_p_dag(f1,2,i1,i),ip
                 J_mu_n_dag(f1,ip,i1,i)=conjg(J_mu_n(f1,ip,i1,i))
              enddo
           enddo
        enddo
     enddo
       res_p=0.0d0
       res_n=0.0d0
       do i1=1,3
          do f1=1,3
             do i=1,4
                do ip=1,6
                   res_p(i,i,ip)=res_p(i,i,ip)+J_mu_p_dag(f1,ip,i1,i)*J_mu_p(f1,ip,i1,i)*xmn**2/(p1(1)*pp1(1))
                   ! write(6,*) res_p(i,i), i, i
                   res_n(i,i,ip)=res_n(i,i,ip)+J_mu_n_dag(f1,ip,i1,i)*J_mu_n(f1,ip,i1,i)*xmn**2/(p1(1)*pp1(1))
                enddo
             enddo
          enddo
       enddo
        
       rl(1, :)=res_p(1,1,:)
       rt(1, :)=res_p(2,2,:)+res_p(3,3,:)
       rl(2, :)=res_n(1,1,:)
       rt(2, :)=res_n(2,2,:)+res_n(3,3,:)
    end subroutine pion_res1b

end module dirac_matrices
