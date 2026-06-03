add_library(gzstream)
target_sources(gzstream
    PRIVATE
    gzstream.C
    )
target_include_directories(gzstream PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(gzstream PRIVATE achilles_options PUBLIC ZLIB::ZLIB)
