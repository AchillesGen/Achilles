add_library(sherpa SHARED
    SherpaMEs.cc
)
target_include_directories(sherpa PUBLIC ${SHERPA_INCLUDE_DIRS})
target_link_libraries(sherpa PRIVATE project_options
                             PUBLIC ${SHERPA_LIBRARIES} spdlog::spdlog)
