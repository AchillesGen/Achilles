module libmap
    use iso_c_binding
    implicit none

    type :: map
        type(c_ptr), private :: dict = c_null_ptr
        contains
            procedure :: init => map_init
            procedure :: exists => map_exists
            procedure :: insert => map_insert
            procedure :: lookup => map_lookup
            procedure :: has => map_contains
            procedure :: empty => map_empty
            procedure :: erase => map_erase
            procedure :: clear => map_clear
            procedure :: delete => map_delete
    end type map

    interface map
        module procedure copy_constructor_map
    end interface

    type :: complex_map
        type(c_ptr), private :: dict
        contains
            procedure :: init => complex_map_init
            procedure :: exists => complex_map_exists
            procedure :: insert => complex_map_insert
            procedure :: lookup => complex_map_lookup
            procedure :: has => complex_map_contains
            procedure :: empty => complex_map_empty
            procedure :: erase => complex_map_erase
            procedure :: clear => complex_map_clear
            procedure :: delete => complex_map_delete
    end type complex_map

    interface complex_map
        module procedure copy_constructor_complex_map
    end interface

    private
    public :: map, complex_map

    interface
        subroutine cmap_init(dict) bind(C, name="map_init")
            use iso_c_binding
            implicit none
            type(c_ptr), intent(out) :: dict
        end subroutine cmap_init

        subroutine cmap_insert(dict, key, val) bind(C, name="map_insert")
            use iso_c_binding
            implicit none
            type(c_ptr), value :: dict
            character(kind=c_char), dimension(*), intent(in) :: key
            real(kind=c_double), intent(in), value :: val
        end subroutine cmap_insert

        function cmap_lookup(dict, key) result(val) bind(C, name="map_lookup")
            use iso_c_binding
            implicit none
            type(c_ptr), value :: dict
            character(kind=c_char), dimension(*), intent(in) :: key
            real(kind=c_double) :: val
        end function cmap_lookup

        function cmap_contains(dict, key) result(val) bind(C, name="map_contains")
            use iso_c_binding
            implicit none
            type(c_ptr), value :: dict
            character(kind=c_char), dimension(*), intent(in) :: key
            logical(kind=c_bool) :: val
        end function cmap_contains        

        function cmap_empty(dict) result(val) bind(C, name="map_empty")
            use iso_c_binding
            implicit none
            type(c_ptr), value :: dict
            logical(kind=c_bool) :: val
        end function cmap_empty

        subroutine cmap_erase(dict, key) bind(C, name="map_erase")
            use iso_c_binding
            implicit none
            type(c_ptr), value :: dict
            character(kind=c_char), dimension(*), intent(in) :: key
        end subroutine cmap_erase

        subroutine cmap_clear(dict) bind(C, name="map_clear")
            use iso_c_binding
            implicit none
            type(c_ptr), value :: dict
        end subroutine cmap_clear

        subroutine cmap_delete(dict) bind(C, name="map_delete")
            use iso_c_binding
            implicit none
            type(c_ptr), value :: dict
        end subroutine cmap_delete

        subroutine ccmap_init(dict) bind(C, name="cmap_init")
            use iso_c_binding
            implicit none
            type(c_ptr), intent(out) :: dict
        end subroutine ccmap_init

        subroutine ccmap_insert(dict, key, val) bind(C, name="cmap_insert")
            use iso_c_binding
            implicit none
            type(c_ptr), value :: dict
            character(kind=c_char), dimension(*), intent(in) :: key
            complex(kind=c_double_complex), intent(in), value :: val
        end subroutine ccmap_insert

        function ccmap_lookup(dict, key) result(val) bind(C, name="cmap_lookup")
            use iso_c_binding
            implicit none
            type(c_ptr), value :: dict
            character(kind=c_char), dimension(*), intent(in) :: key
            complex(kind=c_double_complex), pointer :: val
        end function ccmap_lookup

        function ccmap_contains(dict, key) result(val) bind(C, name="cmap_contains")
            use iso_c_binding
            implicit none
            type(c_ptr), value :: dict
            character(kind=c_char), dimension(*), intent(in) :: key
            logical(kind=c_bool) :: val
        end function ccmap_contains        

        function ccmap_empty(dict) result(val) bind(C, name="cmap_empty")
            use iso_c_binding
            implicit none
            type(c_ptr), value :: dict
            logical(kind=c_bool) :: val
        end function ccmap_empty

        subroutine ccmap_erase(dict, key) bind(C, name="cmap_erase")
            use iso_c_binding
            implicit none
            type(c_ptr), value :: dict
            character(kind=c_char), dimension(*), intent(in) :: key
        end subroutine ccmap_erase

        subroutine ccmap_clear(dict) bind(C, name="cmap_clear")
            use iso_c_binding
            implicit none
            type(c_ptr), value :: dict
        end subroutine ccmap_clear

        subroutine ccmap_delete(dict) bind(C, name="cmap_delete")
            use iso_c_binding
            implicit none
            type(c_ptr), value :: dict
        end subroutine ccmap_delete
    end interface

contains

    subroutine map_init(this)
        class(map), intent(inout) :: this
        call cmap_init(this%dict)
    end subroutine map_init

    function copy_constructor_map(other) result(this)
        type(map) :: this
        type(c_ptr), intent(in), target :: other
        this%dict = other
    end function copy_constructor_map

    function map_exists(this) result(val)
        class(map), intent(in) :: this
        logical :: val
        val = c_associated(this%dict)
    end function map_exists

    subroutine map_insert(this, key, val)
        class(map), intent(inout) :: this
        character(len=*), intent(in) :: key
        real(kind=c_double), intent(in) :: val
        character(len=:), allocatable :: trim_key
        type(c_ptr) :: ckey
        trim_key = trim(key) 
        ckey = f2cstring(trim_key)
        call cmap_insert(this%dict, ckey, val)
    end subroutine map_insert

    function map_lookup(this, key) result(val)
        use libutilities
        class(map), intent(in) :: this
        character(len=*), intent(in) :: key
        character(len=:), allocatable :: trim_key
        type(c_ptr) :: ckey
        real(kind=c_double) :: val
        trim_key = trim(key) 
        ckey = f2cstring(trim_key)
        val = cmap_lookup(this%dict, ckey)
    end function map_lookup

    function map_contains(this, key) result(val)
        class(map), intent(in) :: this
        character(len=*), intent(in) :: key
        character(len=:), allocatable :: trim_key
        type(c_ptr) :: ckey
        logical :: val
        trim_key = trim(key) 
        ckey = f2cstring(trim_key)
        val = cmap_contains(this%dict, ckey)
    end function map_contains

    function map_empty(this) result(val)
        class(map), intent(in) :: this
        logical :: val
        val = cmap_empty(this%dict)
    end function map_empty

    subroutine map_erase(this, key)
        class(map), intent(inout) :: this
        character(len=*), intent(in) :: key
        call cmap_erase(this%dict, key)
    end subroutine map_erase

    subroutine map_clear(this)
        class(map), intent(inout) :: this
        call cmap_clear(this%dict)
    end subroutine map_clear

    subroutine map_delete(this)
        class(map), intent(inout) :: this
        call cmap_delete(this%dict)
    end subroutine map_delete

    subroutine complex_map_init(this)
        class(complex_map), intent(inout) :: this
        call ccmap_init(this%dict)
    end subroutine complex_map_init

    function copy_constructor_complex_map(other) result(this)
        type(complex_map) :: this
        type(c_ptr), intent(in), target :: other
        this%dict = other
    end function copy_constructor_complex_map

    function complex_map_exists(this) result(val)
        class(complex_map), intent(in) :: this
        logical :: val
        val = c_associated(this%dict)
    end function complex_map_exists

    subroutine complex_map_insert(this, key, val)
        class(complex_map), intent(inout) :: this
        character(len=*), intent(in) :: key
        complex(kind=c_double_complex), intent(in) :: val
        character(len=:), allocatable :: trim_key
        type(c_ptr) :: ckey
        trim_key = trim(key) 
        ckey = f2cstring(trim_key)
        call ccmap_insert(this%dict, ckey, val)
    end subroutine complex_map_insert

    function complex_map_lookup(this, key) result(val)
        class(complex_map), intent(in) :: this
        character(len=*), intent(in) :: key
        complex(kind=c_double_complex) :: val
        character(len=:), allocatable :: trim_key
        type(c_ptr) :: ckey
        trim_key = trim(key) 
        ckey = f2cstring(trim_key)
        val = ccmap_lookup(this%dict, ckey)
    end function complex_map_lookup

    function complex_map_contains(this, key) result(val)
        class(complex_map), intent(in) :: this
        character(len=*), intent(in) :: key
        logical :: val
        character(len=:), allocatable :: trim_key
        type(c_ptr) :: ckey
        trim_key = trim(key) 
        ckey = f2cstring(trim_key)
        val = ccmap_contains(this%dict, ckey)
    end function complex_map_contains

    function complex_map_empty(this) result(val)
        class(complex_map), intent(in) :: this
        logical :: val
        val = ccmap_empty(this%dict)
    end function complex_map_empty

    subroutine complex_map_erase(this, key)
        class(complex_map), intent(inout) :: this
        character(len=*), intent(in) :: key
        call ccmap_erase(this%dict, key)
    end subroutine complex_map_erase

    subroutine complex_map_clear(this)
        class(complex_map), intent(inout) :: this
        call ccmap_clear(this%dict)
    end subroutine complex_map_clear

    subroutine complex_map_delete(this)
        class(complex_map), intent(inout) :: this
        call ccmap_delete(this%dict)
    end subroutine complex_map_delete

end module libmap
