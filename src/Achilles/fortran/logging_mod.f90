module liblogging

    use, intrinsic :: iso_c_binding
    implicit none

    private

    type logger_type
        private
        integer :: dummy = 0
    contains
        procedure :: debug => log_debug
        procedure :: info => log_info
        procedure :: error => log_error
        procedure :: warn => log_warn
        procedure :: critical => log_critical
    end type logger_type

    interface
        subroutine debug_c(str) bind(C, name="LogDebug")
            use iso_c_binding
            implicit none
            type(c_ptr), value :: str
        end subroutine

        subroutine info_c(str) bind(C, name="LogInfo")
            use iso_c_binding
            implicit none
            type(c_ptr), value :: str
        end subroutine

        subroutine error_c(str) bind(C, name="LogError")
            use iso_c_binding
            implicit none
            type(c_ptr), value :: str
        end subroutine

        subroutine warn_c(str) bind(C, name="LogWarn")
            use iso_c_binding
            implicit none
            type(c_ptr), value :: str
        end subroutine

        subroutine critical_c(str) bind(C, name="LogCritical")
            use iso_c_binding
            implicit none
            type(c_ptr), value :: str
        end subroutine
    end interface

    type(logger_type) :: logger
    public :: logger

contains
    subroutine log_debug(this, str)
        use libutilities
        implicit none
        class(logger_type) :: this
        character(len=*), intent(in) :: str
        type(c_ptr) :: msg

        msg = f2cstring(str) 
        call debug_c(msg)
        call free(msg)
    end subroutine

    subroutine log_info(this, str)
        use libutilities
        implicit none
        class(logger_type) :: this
        character(len=*), intent(in) :: str
        type(c_ptr) :: msg

        msg = f2cstring(str)
        call info_c(msg)
        call free(msg)
    end subroutine

    subroutine log_error(this, str)
        use libutilities
        implicit none
        class(logger_type) :: this
        character(len=*), intent(in) :: str
        type(c_ptr) :: msg

        msg = f2cstring(str)
        call error_c(msg)
        call free(msg)
    end subroutine

    subroutine log_warn(this, str)
        use libutilities
        implicit none
        class(logger_type) :: this
        character(len=*), intent(in) :: str
        type(c_ptr) :: msg

        msg = f2cstring(str)
        call warn_c(msg)
        call free(msg)
    end subroutine

    subroutine log_critical(this, str)
        use libutilities
        implicit none
        class(logger_type) :: this
        character(len=*), intent(in) :: str
        type(c_ptr) :: msg

        msg = f2cstring(str)
        call critical_c(msg)
        call free(msg)
    end subroutine

end module
