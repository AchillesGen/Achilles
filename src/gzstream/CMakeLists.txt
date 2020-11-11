add_library(gzstream)
target_sources(gzstream
    PRIVATE
    gzstream.C
    )
target_link_libraries(gzstream PRIVATE project_options PUBLIC ZLIB::ZLIB)
