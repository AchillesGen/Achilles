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
