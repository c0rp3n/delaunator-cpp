#pragma once

#if __has_include("delaunator_config.hpp")
#   include "delaunator_config.hpp"
#endif

#ifdef DELAUNATOR_SINGLE_PRECISION
    typedef float dfloat;
#else
    typedef double dfloat;
#endif // DELAUNATOR_SINGLE_PRECISION
