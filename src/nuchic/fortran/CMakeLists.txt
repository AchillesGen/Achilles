enable_language(Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
add_subdirectory(Interface)
add_subdirectory(Tools)
add_subdirectory(Interactions)
add_subdirectory(HardScattering)
