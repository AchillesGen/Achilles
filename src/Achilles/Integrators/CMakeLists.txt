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
list(APPEND achilles_targets AchillesIntegrator)
set(achilles_targets ${achilles_targets} PARENT_SCOPE)
