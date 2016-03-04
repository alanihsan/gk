/*
 * gkgeometry.h
 *
 *  Created on: 2015/06/01
 *      Author: tmakimoto
 */

#ifndef GKGEOMETRY_H_
#define GKGEOMETRY_H_

#include "gkdef.h"

/**
 * @defgroup geometry_tags Geometry Tags
 *
 */

namespace gk {

///*
// * Classification of geometries.
// */
//struct curve_type {
//};
//
//struct surface_type {
//};

/*
 * Classification of traits of a curve.
 */
struct direction_tag {
};

struct vector_tag: public direction_tag {
};

struct curve_tag {
};

struct line_tag: public curve_tag, public direction_tag {
};

struct ray_tag: public line_tag {
};

struct segment_tag: public ray_tag {
};

struct circle_tag: public curve_tag, public direction_tag {
};

struct free_curve_tag: public curve_tag {
};

/*
 * Classification of traits of a surface.
 */

struct surface_tag {
};

struct plane_tag: public surface_tag, public direction_tag {
};

struct sphere_tag: public surface_tag {
};

struct free_surface_tag: public surface_tag {
};

template<typename Category, typename Vector>
struct geometry;

#define GK_GEOMETRY_BASE_TEMPLATE_CLASS(category) \
	template<typename Vector> \
	struct geometry<category, Vector> { \
		typedef category geometry_category; \
		typedef Vector vector_type; \
	};

GK_GEOMETRY_BASE_TEMPLATE_CLASS(direction_tag)
GK_GEOMETRY_BASE_TEMPLATE_CLASS(line_tag)
GK_GEOMETRY_BASE_TEMPLATE_CLASS(ray_tag)
GK_GEOMETRY_BASE_TEMPLATE_CLASS(segment_tag)
GK_GEOMETRY_BASE_TEMPLATE_CLASS(circle_tag)
GK_GEOMETRY_BASE_TEMPLATE_CLASS(free_curve_tag)
GK_GEOMETRY_BASE_TEMPLATE_CLASS(plane_tag)
GK_GEOMETRY_BASE_TEMPLATE_CLASS(sphere_tag)
GK_GEOMETRY_BASE_TEMPLATE_CLASS(free_surface_tag)
#undef GK_GEOMETRY_BASE_TEMPLATE_CLASS

/**
 * @brief Traits of a geometry.
 *
 * @tparam Geometry
 */
template<typename Geometry>
struct geometry_traits {
	typedef typename Geometry::geometry_category geometry_category;
	typedef typename Geometry::vector_type vector_type;
	static const size_t Dimension = Geometry::Dimension;
};

} // namespace gk

#endif /* GKGEOMETRY_H_ */
