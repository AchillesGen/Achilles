module libsystem
    use iso_c_binding

    private
    public :: find_file

    interface
        function find_file_c(filename, head) bind(C, name="FindAchillesFile")
            use iso_c_binding
            implicit none
            type(c_ptr), intent(in), value :: filename, head
            type(c_ptr) :: find_file_c
        end function find_file_c
    end interface

contains

    function find_file(filename, head)
        use libutilities
        character(len=*), intent(in) :: filename, head
        type(c_ptr) :: cfilename, chead, cfile_loc
        character(len=:), allocatable :: find_file

        print*, "here", filename, " ", head
        cfilename = f2cstring(trim(filename))
        chead = f2cstring(trim(head))
        cfile_loc = find_file_c(cfilename, chead)
        find_file = c2fstring(cfile_loc)

    end function find_file

end module libsystem
