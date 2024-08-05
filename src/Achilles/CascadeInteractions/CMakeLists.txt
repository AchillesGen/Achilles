add_library(AchillesCascadeInteractions SHARED
    GeantInteractions.cc
    NasaInteractions.cc
    ConstantInteractions.cc
    MesonBaryonInteractions.cc
    PionAbsorption.cc
    PionInteractions.cc
)
target_link_libraries(AchillesCascadeInteractions PRIVATE project_options project_warnings
                                                  PUBLIC physics)
