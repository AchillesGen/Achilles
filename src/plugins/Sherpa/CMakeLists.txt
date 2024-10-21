add_library(sherpa SHARED
    SherpaInterface.cc
    Channels3.cc
    Channels1.cc
)
if(UFO_VERSION VERSION_GREATER_EQUAL 2.0)
    target_sources(sherpa PRIVATE FormFactors.cc)
    target_compile_definitions(sherpa PUBLIC UFO_2_0)
endif()
target_link_libraries(sherpa PRIVATE project_options
                             PUBLIC SHERPA::All yaml::cpp physics spdlog::spdlog)
if(ENABLE_GZIP)
    target_link_libraries(sherpa PUBLIC ZLIB::ZLIB)
endif()

add_library(SherpaAchillesInput SHARED
    AchillesReader.cc
)
target_link_libraries(SherpaAchillesInput PRIVATE project_options 
                                          PUBLIC SHERPA::All)
