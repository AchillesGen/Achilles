module libspectral_function
    use libinterpolate

    type spectral_function
        private
        type(interp2d) :: pke_p, pke_n
    contains
        procedure :: pke => get_pke
    end type

    interface spectral_function
        module procedure create_spectral_function
    end interface

contains

    function load_spectral_function(fname, mode)
        use liblogging

        implicit none
        character(len=*), intent(in) :: fname, mode
        character(6) :: norm_str
        real*8, parameter :: pi = 4.0*atan(1.0d0)
        type(interp2d) :: load_spectral_function
        real*8, allocatable :: pke(:,:), xe(:), dp(:), p(:)
        real*8 :: hp, he, norm
        integer*4 :: np, ne, i, j

        open(unit=4, file=fname, status='unknown', form='formatted', action='read')
        read(4,*) ne, np
        allocate(p(np), pke(ne,np), dp(np), xe(ne))
        do j=1,np
            read(4,*) p(j)
            read(4,'(4(f6.1,2x,e10.3))')(xe(i),pke(i,j),i=1,ne)
        enddo
        close(4)
        norm=0.0d0
        hp=p(2)-p(1)
        he=xe(2)-xe(1)
        do j=1,np
            dp(j)=sum(pke(:,j))*he
        enddo
        norm=sum(p(:)**2*dp(:))*4.0d0*pi*hp
        write(norm_str,'(f6.4)') norm
        call logger%info("n(k) norm initial for "//mode//"="//norm_str)
        load_spectral_function = interp2d(xe, p, pke, ne, np, 1)
    end function

    function create_spectral_function(fname_pkep, fname_pken)
        implicit none
        character(len=*), intent(in) :: fname_pkep, fname_pken
        type(spectral_function) :: create_spectral_function

        create_spectral_function%pke_p = load_spectral_function(fname_pkep, "proton")
        create_spectral_function%pke_n = load_spectral_function(fname_pken, "neutron")
    end function

    function get_pke(this, mode, e, mom)
        use liblogging

        implicit none
        class(spectral_function), intent(in) :: this
        integer*4, intent(in) :: mode
        real*8, intent(in) :: e, mom
        real*8 :: get_pke
        if(mode.eq.1) then
            get_pke = this%pke_p%call(e, mom)
        else if(mode.eq.2) then
            get_pke = this%pke_n%call(e, mom)
        else
            call logger%critical("Invalid spectral function mode! Use 1 for protons, 2 for neutrons!")
        endif
    end function

endmodule libspectral_function
