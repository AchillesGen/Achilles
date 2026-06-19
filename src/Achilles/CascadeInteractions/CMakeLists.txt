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
target_link_libraries(AchillesCascadeInteractions PRIVATE project_options project_warnings
                                                  PUBLIC physics AchillesIntegrator)
list(APPEND achilles_targets AchillesCascadeInteractions)
set(achilles_targets ${achilles_targets} PARENT_SCOPE)

# The interaction models in this library are only reachable through their static-initializer
# self-registration with the InteractionFactory; no symbol in the library is referenced
# directly. Linkers running with --as-needed (the default on e.g. Ubuntu and Fedora) therefore
# drop the library from the consumer's runtime dependencies, its initializers never run, and
# no cascade interactions are registered. Consumers must link via this variable, which pins
# the library with --no-as-needed and restores the previous linker state afterwards.
if(UNIX AND NOT APPLE)
    set(AchillesCascadeInteractionsLink
        "-Wl,--push-state,--no-as-needed" AchillesCascadeInteractions "-Wl,--pop-state"
        PARENT_SCOPE)
else()
    set(AchillesCascadeInteractionsLink AchillesCascadeInteractions PARENT_SCOPE)
endif()
