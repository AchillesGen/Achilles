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

# gzstream/gzstream.h is included from the public EventWriter.hh, so it must ship
# with the install for downstream consumers.
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/gzstream
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
