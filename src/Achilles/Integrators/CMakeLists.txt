add_library(AchillesIntegrator SHARED
    AdaptiveQuadrature.cc
    Adzint.cc
    BaseIntegrator.cc
    DoubleExponential.cc
    GaussKronrod.cc
    NewtonCotes.cc
)
target_link_libraries(AchillesIntegrator PRIVATE project_warnings
                                         PUBLIC utilities)
get_directory_property(SUBDIR_TARGETS DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DEFINITON achilles_targets)
set(achilles_targets ${achilles_targets} AchillesIntegrator PARENT_SCOPE)
