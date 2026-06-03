add_library(AchillesWriters SHARED EventWriter.cc)
target_link_libraries(AchillesWriters PRIVATE achilles_options achilles_warnings
                                      PUBLIC fmt::fmt spdlog::spdlog yaml-cpp::yaml-cpp
                                             utilities physics nuclear_models)

if(ACHILLES_ENABLE_HEPMC3)
    target_sources(AchillesWriters PRIVATE
        HepMC3EventWriter.cc
        NuHepMCWriter.cc
    )
    target_link_libraries(AchillesWriters PUBLIC HepMC3::HepMC3 NuHepMC::CPPUtils)
    # Treat HepMC3 / NuHepMC headers as system headers so their warnings don't trip -Werror.
    target_include_directories(AchillesWriters SYSTEM PRIVATE
        $<TARGET_PROPERTY:HepMC3::HepMC3,INTERFACE_INCLUDE_DIRECTORIES>
        $<TARGET_PROPERTY:NuHepMC::CPPUtils,INTERFACE_INCLUDE_DIRECTORIES>)
endif()

list(APPEND achilles_targets AchillesWriters)
set(achilles_targets ${achilles_targets} PARENT_SCOPE)
