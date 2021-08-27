add_library(nuchic_integrators SHARED
    AdaptiveIntegrator.cc
    BaseIntegrator.cc
    DoubleExponential.cc
    GaussKronrod.cc
    NewtonCotes.cc
)
target_link_libraries(nuchic_integrators PRIVATE project_options project_warnings
                                         PUBLIC spdlog::spdlog)
