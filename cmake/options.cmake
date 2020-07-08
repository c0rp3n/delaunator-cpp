option(WITH_TESTS
    "Choose if delaunator tests should be built" TRUE)
add_feature_info("Unit tests" WITH_TESTS "Delaunator unit tests")

option(DELAUNATOR_HEADER_ONLY "Whether to the library should be header only." OFF)
option(DELAUNATOR_SINGLE_PRECISION "Whether to use single precision floating point values." OFF)
