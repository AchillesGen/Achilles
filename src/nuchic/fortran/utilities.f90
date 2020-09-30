module libutilities
    use, intrinsic :: iso_c_binding
    implicit none

    private

    public :: constants, c2fstring, f2cstring, malloc, free

    type constants_type
        double precision :: c, hbarc, hbarc2, mp, mn, mqe
    end type constants_type

    interface
        function strlen(str) result(iszie) bind(C, name="strlen")
            use iso_c_binding
            type(c_ptr), value :: str
            integer(c_int) :: isize
        end function strlen

        function malloc(isize) bind(C, name="malloc")
            use iso_c_binding
            type(c_ptr) :: malloc
            integer(c_int), value, intent(in) :: isize
        end function malloc

        subroutine free(ptr) bind(C, name="free")
            use iso_c_binding
            type(c_ptr), value, intent(in) :: ptr
        end subroutine free
    end interface

    type(constants_type) :: constants

contains

    function c2fstring(cstr) result(fstr)
        character (len=:), allocatable :: fstr

        type(c_ptr) :: cstr
        integer(c_int) :: n

        if(c_associated(cstr)) then
            n = strlen(cstr)

            block
                ! Convert the C string to a Fortran string
                character(kind=c_char, len=n+1), pointer :: s
                call c_f_pointer(cptr=cstr,fptr=s)
                fstr = s(1:n)
                nullify(s)
            end block
        else
            fstr = ''
        end if
    end function

    function f2cstring(fstr) result(cstr)
        character(len=*), intent(in) :: fstr
        type(c_ptr) :: cstr
        character(len=1,kind=c_char), pointer :: cptr(:)
        integer :: i, length

        length = len(fstr)
        if(length <= 0) then
            cstr = c_null_ptr
        else
            cstr = malloc(length+1)
            if(c_associated(cstr)) then
                call c_f_pointer(cstr, cptr,[length+1])
                forall(i=1:length)
                    cptr(i) = fstr(i:i)
                end forall
                cptr(length+1) = c_null_char
            endif
        endif
    end function f2cstring

end module
