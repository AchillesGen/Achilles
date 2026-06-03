add_library(AchillesCascadeInteractions SHARED
    GeantInteractions.cc
    NasaInteractions.cc
    ConstantInteractions.cc
    MesonBaryonInteractions.cc
    OsetMesonBaryonInteractions.cc
    PionAbsorption.cc
    PionInteractions.cc
    DeltaInteractions.cc
    NucleonNucleon.cc
)
target_link_libraries(AchillesCascadeInteractions PRIVATE achilles_options achilles_warnings HighFive
                                                  PUBLIC physics AchillesIntegrator)
list(APPEND achilles_targets AchillesCascadeInteractions)
set(achilles_targets ${achilles_targets} PARENT_SCOPE)
