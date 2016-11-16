/*
 * gkgeometry.h
 *
 *  Created on: 2015/06/01
 *      Author: tmakimoto
 */

#ifndef GKGEOMETRY_H_
#define GKGEOMETRY_H_

#include "gkdef.h"
#include <Eigen/Core>

/**
 * @defgroup geometry_tags Geometry Tags
 *
 */

namespace gk {

/*
 * Classification of traits of a curve.
 */

/**
 * @brief Direction category.
 */
struct direction_tag {
};

/**
 * @brief Curve category.
 */
struct curve_tag {
};

/**
 * @brief Line category.
 */
struct line_tag: public curve_tag, public direction_tag {
};

/**
 * @brief Ray category.
 */
struct ray_tag: public line_tag {
};

/**
 * @brief Segment category.
 */
struct segment_tag: public ray_tag {
};

/**
 * @brief Circle category.
 */
struct circle_tag: public curve_tag, public direction_tag {
};

/**
 * @brief Free curve category.
 */
struct free_curve_tag: public curve_tag {
};

/*
 * Classification of traits of a surface.
 */

/**
 * @brief Surface category.
 */
struct surface_tag {
};

/**
 * @brief Plane category.
 */
struct plane_tag: public surface_tag, public direction_tag {
};

/**
 * @brief Sphere category.
 */
struct sphere_tag: public surface_tag {
};

template<typename Category, typename T, std::size_t Dimension>
struct geometry;

#define GK_GEOMETRY_BASE_TEMPLATE_CLASS(Category) \
	template<typename T, std::size_t Dimension> \
	struct geometry<Category, T, Dimension> { \
		typedef Category category; \
		typedef T value_type; \
		typedef Eigen::Matrix<T, 1, Dimension, Eigen::RowMajor> vector_type; \
};

GK_GEOMETRY_BASE_TEMPLATE_CLASS(direction_tag)
GK_GEOMETRY_BASE_TEMPLATE_CLASS(line_tag)
GK_GEOMETRY_BASE_TEMPLATE_CLASS(ray_tag)
GK_GEOMETRY_BASE_TEMPLATE_CLASS(segment_tag)
GK_GEOMETRY_BASE_TEMPLATE_CLASS(circle_tag)
GK_GEOMETRY_BASE_TEMPLATE_CLASS(free_curve_tag)
GK_GEOMETRY_BASE_TEMPLATE_CLASS(plane_tag)
GK_GEOMETRY_BASE_TEMPLATE_CLASS(sphere_tag)

#undef GK_GEOMETRY_BASE_TEMPLATE_CLASS

/**
 * @brief Traits of a geometry.
 *
 * @tparam Geometry
 */
template<typename Geometry>
struct geometry_traits {
	typedef typename Geometry::category category;
	typedef typename Geometry::vector_type vector_type;
	static const size_t Dimension = Geometry::Dimension;
};

} // namespace gk

#endif /* GKGEOMETRY_H_ */
