enable_language(Fortran)

add_library(fortran_interface SHARED
    utilities.f90
    vectors_mod.f90
    particle_mod.f90
    particle_info_mod.f90
    logging_mod.f90
    interaction_mod.f90
    interpolation_mod.f90
    # hardscattering_mod.f90
    event_mod.f90
    nucleus_mod.f90
    nuclear_model.f90
    nuclear_model_factory.f90
    nuclear_model_interface.f90
    process_info_mod.f90
    OneBody/pke_read_xsec.f90
    OneBody/contract.f90
    OneBody/currents.f90
    OneBody/nform.f90
    OneBody/mathtool.f90
    # )
    # target_link_libraries(fortran_interface_f PRIVATE project_options
    #                                           PUBLIC physics)
    # 
    # add_library(fortran_interface SHARED
    FConstants.cc
    FVectors.cc
    FParticle.cc
    FParticleInfo.cc
    FLogging.cc
    FInteractions.cc
    FInterpolation.cc
    # FHardScattering.cc
    FEvent.cc 
    FNucleus.cc
    FProcessInfo.cc
    FNuclearModel.cc
)
target_link_libraries(fortran_interface PRIVATE project_options #$<BUILD_INTERFACE:project_warnings>
                                        PUBLIC physics nuclear_models)

add_library(AchillesPlugin_nuclear_test nuclear_model_test.f90)
target_link_libraries(AchillesPlugin_nuclear_test PRIVATE project_options
                                                  PUBLIC fortran_interface)

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
