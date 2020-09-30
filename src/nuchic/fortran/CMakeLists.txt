enable_language(Fortran)

add_executable(constant_test 
    utilities.f90
    test_constants.f90
    FConstants.cc
    FVectors.cc
    FParticle.cc
    FParticleInfo.cc
    FLogging.cc
    vectors_mod.f90
    particle_mod.f90
    particle_info_mod.f90
    logging_mod.f90
    )
target_link_libraries(constant_test PRIVATE project_options
                                    PUBLIC physics)
