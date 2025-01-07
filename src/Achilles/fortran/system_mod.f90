module libsystem

    use, intrinsic :: iso_c_binding
    implicit none

    private

    interface
        function find_file_c(file, head) bind(C, name="FindFile")
            use iso_c_binding
            implicit none
            type(c_ptr), value :: file, head
            type(c_ptr) :: find_file_c
        end function
    end interface

    public :: find_file

contains
    function find_file(file, head)
        use libutilities
        implicit none
        type(c_ptr) :: cfile, chead
        character(len=*), intent(in) :: file, head
        character(len=:), allocatable :: find_file

        cfile = f2cstring(trim(file))
        chead = f2cstring(trim(head))
        find_file = c2fstring(find_file_c(cfile, chead))
        call free(cfile)
        call free(chead)
    end function
end module
