# Third-party; linked statically and privately so it is neither installed nor exported
add_library(gzstream STATIC)
target_sources(gzstream
    PRIVATE
    gzstream.C
    )
set_target_properties(gzstream PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(gzstream PRIVATE achilles_options PUBLIC ZLIB::ZLIB)
