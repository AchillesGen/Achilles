add_library(fortran_tools SHARED
    spectral_function.f90
    lorentz.f90
    nform.f90
)
target_link_libraries(fortran_tools PRIVATE project_options
                                    PUBLIC fortran_interface)
