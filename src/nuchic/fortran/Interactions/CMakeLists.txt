add_library(fortran_interactions SHARED
    interaction_mod.f90
    FInteractionComponent.cc
)
target_link_libraries(fortran_interactions PRIVATE project_options #$<BUILD_INTERFACE:project_warnings>
                                           PUBLIC fortran_tools fortran_interface)
