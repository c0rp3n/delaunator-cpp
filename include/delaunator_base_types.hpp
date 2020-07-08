#pragma once

#include "delaunator_config.hpp"

#ifdef DELAUNATOR_SINGLE_PRECISION
    typedef float dfloat;
#else
    typedef double dfloat;
#endif // DELAUNATOR_SINGLE_PRECISION
