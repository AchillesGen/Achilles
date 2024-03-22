add_library(sherpa SHARED
    SherpaInterface.cc
    Channels3.cc
    Channels1.cc
)
if(UFO_VERSION VERSION_GREATER_EQUAL 2.0)
    target_sources(sherpa FormFactors.cc)
endif()
target_include_directories(sherpa PUBLIC ${SHERPA_INCLUDE_DIRS})
target_link_libraries(sherpa PRIVATE project_options
                      PUBLIC ${SHERPA_LIBRARIES} physics spdlog::spdlog yaml::cpp)
if(ENABLE_GZIP)
    target_link_libraries(sherpa PUBLIC ZLIB::ZLIB)
endif()
