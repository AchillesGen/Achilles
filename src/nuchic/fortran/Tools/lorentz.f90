module liblorentz

    private
        complex*16, parameter :: czero = (0.0d0, 0.0d0)
        complex*16, parameter :: cone  = (1.0d0, 0.0d0)
        complex*16, parameter :: ci    = (0.0d0, 1.0d0)

    type lorentz_type
            complex*16 :: sig(3, 2, 2), id(2, 2), id4(4, 4), up(2), down(2)
            complex*16 :: gamma_mu(4, 4, 5), g_munu(4, 4), sigma_munu(4, 4, 4, 4)
        contains
            procedure :: initialize => lorentz_init
            procedure :: spinor => lorentz_spinor
    end type lorentz_type

    type(lorentz_type) :: lorentz
    public :: lorentz
    
contains

    subroutine lorentz_init(this)
        implicit none
        logical, save :: initialized = .false.
        class(lorentz_type) :: this
        integer*4 :: i, j

        if(initialized) then
            return
        endif

        ! I_2 initialization
        this%id(:,:)=czero
        this%id(1,1)=cone
        this%id(2,2)=cone

        ! Pauli matrix initialization
        this%sig(:,:,:)=czero
        this%sig(1,1,2)=cone
        this%sig(1,2,1)=cone
        this%sig(2,1,2)=-ci
        this%sig(2,2,1)=ci
        this%sig(3,1,1)=cone
        this%sig(3,2,2)=-cone


        ! I_4 initialization
        this%id4=czero    
        this%id4(1:2,1:2)=this%id
        this%id4(3:4,3:4)=this%id

        ! Gamma matrix initialization
        this%gamma_mu=czero
        this%gamma_mu(1:2,1:2,1)=this%id
        this%gamma_mu(3:4,3:4,1)=-this%id
        do i=2,4
          this%gamma_mu(1:2,3:4,i)=this%sig(i-1,:,:)
          this%gamma_mu(3:4,1:2,i)=-this%sig(i-1,:,:)
        enddo
        this%gamma_mu(1:2,3:4,5)=this%id
        this%gamma_mu(3:4,1:2,5)=this%id

        ! Metric initialization
        this%g_munu=czero
        this%g_munu(1,1)=cone
        this%g_munu(2,2)=-cone
        this%g_munu(3,3)=-cone
        this%g_munu(4,4)=-cone

        ! Sigma_\mu\nu initialization
        do i=1,4
           do j=1,4
              this%sigma_munu(:,:,i,j)=ci*0.5d0*(matmul(this%gamma_mu(:,:,i),this%gamma_mu(:,:,j)) &
                   &     -matmul(this%gamma_mu(:,:,j),this%gamma_mu(:,:,i)))
           enddo
        enddo
        
        ! Spinor initialization
        this%up(1)=cone;this%up(2)=czero
        this%down(1)=czero;this%down(2)=cone

        initialized = .true.
    end subroutine

    subroutine lorentz_spinor(this, mom, mass, spinorU, spinorubar) 
        implicit none
        class(lorentz_type), intent(in) :: this
        real*8, intent(in) :: mom(4), mass
        complex*16, intent(out) :: SpinorU(2, 4), SpinorUbar(2, 4)
        integer*4 :: i
        complex*16 :: sigp(2, 2)
        real*8 :: norm
        sigp = czero
        ! Initialize spinor
        SpinorU=czero
        SpinorUbar=czero
        ! Initialize normalization factors
        norm = sqrt((mom(1)+mass)/(2.0*mom(1)))
        ! Define simga*p
        do i = 1,3
            sigp = sigp + this%sig(i,:,:)*mom(i+1)
        enddo
        ! Build spinor
        SpinorU(1, 1:2) = this%up(:)
        SpinorU(1, 3:4) = matmul(sigp(:, :), this%up(:))/(mom(1)+mass)
        SpinorU(2, 1:2) = this%down(:)
        SpinorU(2, 3:4) = matmul(sigp(:, :), this%down(:))/(mom(1)+mass)
        SpinorU(:, :) = norm*SpinorU(:, :)
        
        SpinorUbar(1, 1:2) = this%up(:)
        SpinorUbar(1, 3:4) = -matmul(this%up(:), sigp(:, :))/(mom(1)+mass)
        SpinorUbar(2, 1:2) = this%down(:)
        SpinorUbar(2, 3:4) = -matmul(this%down(:), sigp(:, :))/(mom(1)+mass)
        SpinorUbar(:, :) = norm*SpinorUbar(:, :)
    end subroutine lorentz_spinor

end module liblorentz
