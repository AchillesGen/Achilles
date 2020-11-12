program test_yaml

   use yaml_version, only: yaml_commit_id=>git_commit_id, &
                           yaml_branch_name=>git_branch_name
   use yaml_types
   use yaml
   use, intrinsic :: iso_fortran_env

   character(error_length) :: error
   character(256) :: path
   class (type_node),pointer :: root

   write(*,*)
   write(*,*) 'YAML version:   ',yaml_commit_id,' (',yaml_branch_name,' branch)'
   write(*,*)

   call get_command_argument(1, path)
   if (path=='') then
      write (*,*) 'ERROR: path to YAML file not provided.'
      stop 2
   end if
   root => parse(path,unit=100,error=error)
   if (error/='') then
      write (*,*) 'PARSE ERROR: '//trim(error)
      stop 1
   end if
   call root%dump(unit=output_unit,indent=0)

end program
