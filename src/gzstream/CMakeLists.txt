add_library(gzstream)
target_sources(gzstream
    PRIVATE
    gzstream.C
    )
target_link_libraries(gzstream PRIVATE achilles_options PUBLIC ZLIB::ZLIB)
