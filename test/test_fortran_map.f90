module test_map_interface
    use libmap

    private
    public :: test_map
    type(map) :: mapping
    type(complex_map) :: complex_mapping

contains

    subroutine test_map(test_suite)
        use unit_test
        type(test_suite_type) :: test_suite
        type(map) :: mapping
        character(len=100) :: key

        call test_case_create('Double Map interface test', test_suite)
        call assert_true(test_map_init(mapping), suite=test_suite)
        call assert_true(test_map_insert_lookup(mapping, 'key1', 1d0), suite=test_suite)

        call test_case_create('Complex Map interface test', test_suite)
        call assert_true(test_cmap_init(complex_mapping), suite=test_suite)
        call assert_true(test_cmap_insert_lookup(complex_mapping, 'key1', (1d0, 1d0)), suite=test_suite)
    end subroutine test_map

    function test_map_init(mapping) result(res)
        type(map) :: mapping
        logical :: res

        res = mapping%exists()
        if(res .eqv. .true.) then
            res = .false.
            return 
        end if
        call mapping%init()
        res = mapping%exists()
    end function test_map_init

    function test_map_insert_lookup(mapping, key, val) result(res)
        type(map) :: mapping
        character(len=*) :: key
        real(8) :: val
        logical :: res

        res = mapping%exists()
        if(res .eqv. .false.) then
            call mapping%init()
        end if

        call mapping%insert(key, val)
        res = mapping%lookup(key) == val
    end function test_map_insert_lookup

    function test_cmap_init(mapping) result(res)
        type(complex_map) :: mapping
        logical :: res

        res = mapping%exists()
        if(res .eqv. .true.) then
            res = .false.
            return 
        end if
        call mapping%init()
        res = mapping%exists()
    end function test_cmap_init

    function test_cmap_insert_lookup(mapping, key, val) result(res)
        type(complex_map) :: mapping
        character(len=*) :: key
        complex(8) :: val
        logical :: res

        res = mapping%exists()
        if(res .eqv. .false.) then
            call mapping%init()
        end if

        call mapping%insert(key, val)
        print*, val, mapping%lookup(key)
        res = mapping%lookup(key) == val
    end function test_cmap_insert_lookup

end module test_map_interface
