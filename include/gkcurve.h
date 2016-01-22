/*
 * gkcurve.h
 *
 *  Created on: 2015/09/25
 *      Author: makitaku
 */

#ifndef GKCURVE_H_
#define GKCURVE_H_

#include "gkgeometry.h"
#include "gkvector.h"
#include "gkaabb.h"

namespace gk {

template<typename CurveCategory, typename Vector, typename Parameter>
struct curve;

#define GK_CURVE_BASE_TEMPLATE_CLASS(category) \
		template<typename Vector, typename Parameter> \
		struct curve<category, Vector, Parameter> : public geometry<category, Vector, Parameter> { \
		};

GK_CURVE_BASE_TEMPLATE_CLASS(curve_tag)
GK_CURVE_BASE_TEMPLATE_CLASS(line_tag)
GK_CURVE_BASE_TEMPLATE_CLASS(ray_tag)
GK_CURVE_BASE_TEMPLATE_CLASS(segment_tag)
GK_CURVE_BASE_TEMPLATE_CLASS(circle_tag)
GK_CURVE_BASE_TEMPLATE_CLASS(free_curve_tag)
#undef GK_CURVE_BASE_TEMPLATE_CLASS

template<typename Curve>
struct curve_traits: public geometry_traits<Curve> {
	typedef geometry_traits<Curve> base;
	typedef typename base::vector_type vector_type;
	typedef typename base::parameter parameter;
	typedef typename vector_traits<vector_type>::value_type value_type;
	typedef value_type distance_type;
	typedef value_type length_type;
	typedef std::pair<vector_type, vector_type> segment_type;
	typedef aabb<vector_type> boundary_type;
};

/**
 * @brief Returns true if a curve is closed, otherwise false.
 * @param curve
 * @return
 */
template<typename Curve>
bool is_closed(const Curve& curve);

/**
// * @brief Returns a domain of a curve.
// * @param curve
// * @return
// */
//template<typename Curve>
//std::pair<typename curve_traits<Curve>::parameter,
//		typename curve_traits<Curve>::parameter> domain(const Curve& curve);

/**
 * @brief
 * @param x
 * @return
 */
template<typename Curve>
typename curve_traits<Curve>::distance_type length(const Curve& x);

template<typename Curve, typename T>
typename divides_result<Curve, T>::value_type derivative(const Curve& x);

}  // namespace gk

#endif /* GKCURVE_H_ */
