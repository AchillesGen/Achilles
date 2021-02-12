add_library(fortran_hardscattering SHARED
    contract.f90
    hardscattering_mod.f90
    hardscattering_calculator.f90
    onebody.f90
    pion.f90
    amp_dcc_sl.f
    FHardScattering.cc
)
target_link_libraries(fortran_hardscattering PRIVATE project_options #$<BUILD_INTERFACE:project_warnings>
                                             PUBLIC fortran_tools fortran_interface)
