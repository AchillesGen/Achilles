add_library(gzstream)
target_sources(gzstream
    PRIVATE
    gzstream.C
    )
target_link_libraries(gzstream PRIVATE achilles_options PUBLIC ZLIB::ZLIB)

# Included by the public EventWriter.hh
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/gzstream
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
