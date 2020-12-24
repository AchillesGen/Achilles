enable_language(Fortran)

add_library(fortran_interface_f SHARED
    utilities.f90
    vectors_mod.f90
    particle_mod.f90
    particle_info_mod.f90
    logging_mod.f90
    interaction_mod.f90
    interpolation_mod.f90
    hardscattering_mod.f90
    OneBody/pke_read_xsec.f90
    OneBody/contract.f90
    OneBody/currents.f90
    OneBody/nform.f90
    OneBody/mathtool.f90
)
target_link_libraries(fortran_interface_f PRIVATE project_options
                                          PUBLIC physics)

add_library(fortran_interface SHARED
    FConstants.cc
    FVectors.cc
    FParticle.cc
    FParticleInfo.cc
    FLogging.cc
    FInteractions.cc
    FInterpolation.cc
    FHardScattering.cc
)
target_link_libraries(fortran_interface PRIVATE project_options $<BUILD_INTERFACE:project_warnings>
                                        PUBLIC fortran_interface_f)


# add_executable(constant_test 
#     test_constants.f90
# )
# target_link_libraries(constant_test PRIVATE project_options
#                                     PUBLIC fortran_interface
# )
# 
# add_executable(interaction_test
#     test_interaction.cc
# )
# target_link_libraries(interaction_test PRIVATE project_options
#                                        PUBLIC fortran_interface
# )
# 
# add_executable(fqe_test test_fqe.cc)
# target_link_libraries(fqe_test PRIVATE project_options
#                                PUBLIC fortran_interface
# )
