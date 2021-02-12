add_library(fortran_interface SHARED
    utilities.f90
    vectors_mod.f90
    particle_mod.f90
    particle_info_mod.f90
    logging_mod.f90
    interpolation_mod.f90
    FConstants.cc
    FVectors.cc
    FParticle.cc
    FParticleInfo.cc
    FLogging.cc
    FInterpolation.cc
)
target_link_libraries(fortran_interface PRIVATE project_options #$<BUILD_INTERFACE:project_warnings>
                                        PUBLIC physics)
