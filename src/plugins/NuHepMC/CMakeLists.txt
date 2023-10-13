CPMAddPackage(
  NAME NuHepMC_CPPUtils
  GIT_TAG main
  GIT_REPOSITORY "https://github.com/NuHepMC/cpputils.git"
  OPTIONS "BUILTIN_HEPMC3 OFF"
)

add_library(nuhepmc SHARED NuHepMCWriter.cc)
target_link_libraries(nuhepmc PRIVATE project_options
                               PUBLIC NuHepMC::CPPUtils fmt::fmt spdlog::spdlog yaml::cpp utilities physics)
