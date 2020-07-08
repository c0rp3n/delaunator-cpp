#pragma once

#include "delaunator_config.hpp"

#ifdef DELAUNATOR_HEADER_ONLY
#   define DELAUNATOR_CONFIGURABLE
#   if __has_include(<concepts>)
//#       define DELAUNATOR_HAS_CONCEPTS
#   endif // __has_include(<concepts>)
#endif // DELAUNATOR_SINGLE_HEADER

#ifdef DELAUNATOR_CONFIGURABLE
#   define DELAUNATOR_MTEMPLATE template<class Config>
#   ifndef DELAUNATOR_HAS_CONCEPTS
#       define DELAUNATOR_TEMPLATE template<class Config = DefaultPointConfig>
#   else
#       define DELAUNATOR_TEMPLATE template <class Config = DefaultPointConfig>\
requires delaunator::PointConfigConcept<Config>

#       ifndef __cpp_lib_concepts
#           define __cpp_lib_concepts
#       endif // __cpp_lib_concepts

#       include <concepts>

#       include "delaunator_base_types.hpp"

        namespace delaunator
        {
            template<class PointConfig>
            concept PointConfigConcept = requires(PointConfig::point_type p0, PointConfig::point_type p1)
            {
                PointConfig::point_type;
                PointConfig::point_type(0.0, 0.0);
                PointConfig::get_x(p0);
                PointConfig::get_y(p0);
                PointConfig::get_magnitude2(p0);
                PointConfig::get_determinant(p0, p1);
                PointConfig::get_vector(p0, p1);
                PointConfig::get_dist2(p0, p1);
                PointConfig::get_equal(p0, p1, 1.0);
            };

            template <typename T, typename U>
            concept Container = requires(T a)
            {
                std::data(a);
                std::size(a);
                std::same_as<typename T::value_type, U>;
            };
        }
#   endif // DELAUNATOR_HAS_CONCEPTS
#else
#   define DELAUNATOR_TEMPLATE
#   define DELAUNATOR_MTEMPLATE
#endif // DELAUNATOR_CONFIGURABLE
