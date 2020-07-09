#pragma once

#if __has_include("delaunator_config.hpp")
#   include "delaunator_config.hpp"
#endif
#include "delaunator_base_types.hpp"
#include "delaunator_concepts.hpp"

#ifdef DELAUNATOR_HEADER_ONLY
#   define INLINE inline
#else
#   define INLINE
#endif // DELAUNATOR_HEADER_ONLY

#include <limits>
#include <vector>
#include <ostream>

#ifndef DELAUNATOR_CONFIGURABLE
#   define DELAUNATOR_CLASS Delaunator
#else
#   define DELAUNATOR_CLASS Delaunator<Config>
#endif

namespace delaunator {

constexpr std::size_t INVALID_INDEX =
    (std::numeric_limits<std::size_t>::max)();

class Point
{
public:
    Point(dfloat x, dfloat y) : m_x(x), m_y(y)
    {}
    Point() : m_x(0), m_y(0)
    {}

    dfloat x() const
    { return m_x; }

    dfloat y() const
    { return m_y; }

    dfloat magnitude2() const
    { return m_x * m_x + m_y * m_y; }

    static dfloat determinant(const Point& p1, const Point& p2)
    {
        return p1.m_x * p2.m_y - p1.m_y * p2.m_x;
    }

    static Point vector(const Point& p1, const Point& p2)
    {
        return Point(p2.m_x - p1.m_x, p2.m_y - p1.m_y);
    }

    static dfloat dist2(const Point& p1, const Point& p2)
    {
        Point vec = vector(p1, p2);
        return vec.m_x * vec.m_x + vec.m_y * vec.m_y;
    }

    static bool equal(const Point& p1, const Point& p2, dfloat span)
    {
        dfloat dist = dist2(p1, p2) / span;

        // ABELL - This number should be examined to figure how how
        // it correlates with the breakdown of calculating determinants.
        return dist < 1e-20;
    }

private:
    dfloat m_x;
    dfloat m_y;
};

class DefaultPointConfig
{
public:
    using point_type = Point;

    static inline dfloat get_x(const point_type& p);
    static inline dfloat get_y(const point_type& p);
    static inline dfloat get_magnitude2(const point_type& p);
    static inline dfloat get_determinant(const point_type& p0, const point_type& p1);
    static inline point_type get_vector(const point_type& p0, const point_type& p1);
    static inline dfloat get_dist2(const point_type& p0, const point_type& p1);
    static inline bool get_equal(const point_type& p0, const point_type& p1, dfloat span);
};

#ifndef DELAUNATOR_CONFIGURABLE
    typedef DefaultPointConfig Config;
#endif

inline std::ostream& operator<<(std::ostream& out, const Point& p)
{
    out << p.x() << "/" << p.y();
    return out;
}

DELAUNATOR_TEMPLATE
class Points
{
public:
    using point_type = typename Config::point_type;
    using const_iterator = point_type const *;

    Points(const std::vector<dfloat>& coords)
        : m_data(reinterpret_cast<const point_type*>(std::data(coords))),
          m_size(std::size(coords) >> 1)
    {}

#   ifndef DELAUNATOR_HAS_CONCEPTS
    template <class Array>
#   else
    template <class Array> requires Container<Array, point_type>
#   endif
    Points(const Array& coords)
        : m_data(std::data(coords)),
          m_size(std::size(coords))
        { 
            static_assert(std::is_same<typename Array::value_type, point_type>::value,
                          "input point type does not match the configured "
                          "point type.")
        }

    const point_type& operator[](size_t offset)
        { return *(m_data + offset); };

    Points::const_iterator begin() const
        { return m_data; }
    Points::const_iterator end() const
        { return m_data + m_size; }
    size_t size() const
        { return m_size; }

private:
    const point_type* m_data;
    const size_t m_size;
};

DELAUNATOR_TEMPLATE
class Delaunator {

public:
    using point_type = typename Config::point_type;

#ifdef DELAUNATOR_CONFIGURABLE
    Points<Config> m_points;
#else
    Points m_points;
#endif

    // 'triangles' stores the indices to the 'X's of the input 'coords'.
    std::vector<std::size_t> triangles;

    // 'halfedges' store indices into 'triangles'.  If halfedges[X] = Y,
    // It says that there's an edge from X to Y where a) X and Y are
    // both indices into triangles and b) X and Y are indices into different
    // triangles in the array.  This allows you to get from a triangle to
    // its adjacent triangle.  If the a triangle edge has no adjacent triangle,
    // its half edge will be INVALID_INDEX.
    std::vector<std::size_t> halfedges;

    std::vector<std::size_t> hull_prev;
    std::vector<std::size_t> hull_next;

    // This contains indexes into the triangles array.
    std::vector<std::size_t> hull_tri;
    std::size_t hull_start;

#   ifndef DELAUNATOR_HAS_CONCEPTS
    template<class Array>
#   else
    template<class Array> requires Container<Array, point_type>
#   endif
    inline Delaunator(Array const& in_coords);
    INLINE Delaunator(std::vector<dfloat> const& in_coords);
    INLINE dfloat get_hull_area();
    INLINE dfloat get_triangle_area();

private:
    std::vector<std::size_t> m_hash;
    point_type m_center;
    std::size_t m_hash_size;
    std::vector<std::size_t> m_edge_stack;

    INLINE void triangulate();
    INLINE std::size_t legalize(std::size_t a);
    INLINE std::size_t hash_key(const point_type& p) const;
    INLINE std::size_t add_triangle(
        std::size_t i0,
        std::size_t i1,
        std::size_t i2,
        std::size_t a,
        std::size_t b,
        std::size_t c);
    INLINE void link(std::size_t a, std::size_t b);
};

} //namespace delaunator

dfloat delaunator::DefaultPointConfig::get_x(const point_type& p)
{
    return p.x();
}

dfloat delaunator::DefaultPointConfig::get_y(const point_type& p)
{
    return p.y();
}

dfloat delaunator::DefaultPointConfig::get_magnitude2(const point_type& p)
{
    return p.magnitude2();
}

dfloat delaunator::DefaultPointConfig::get_determinant(const point_type& p0, const point_type& p1){
    return point_type::determinant(p0, p1);
}

delaunator::DefaultPointConfig::point_type
    delaunator::DefaultPointConfig::get_vector(const point_type& p0, const point_type& p1) {
    return point_type::vector(p0, p1);
}

dfloat delaunator::DefaultPointConfig::get_dist2(const point_type& p0, const point_type& p1) {
    return point_type::dist2(p0, p1);
}

bool delaunator::DefaultPointConfig::get_equal(const point_type& p0, const point_type& p1, dfloat span) {
    return point_type::equal(p0, p1, span);
}

DELAUNATOR_MTEMPLATE
#ifndef DELAUNATOR_HAS_CONCEPTS
template<class Array>
#else
template<class Array> requires Container<Array, point_type>
#endif
delaunator::DELAUNATOR_CLASS::Delaunator(Array const& in_coords)
    : m_points(in_coords)
{
    this->triangulate();
}

#ifdef DELAUNATOR_HEADER_ONLY
#   include "delaunator.cpp"
#endif

#undef INLINE
