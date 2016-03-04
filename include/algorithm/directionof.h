/*
 * directionof.h
 *
 *  Created on: 2016/01/24
 *      Author: tmakimoto
 */

#ifndef INCLUDE_ALGORITHM_DIRECTIONOF_H_
#define INCLUDE_ALGORITHM_DIRECTIONOF_H_

#include "../gkvector.h"

namespace gk {

template<typename Geometry, typename Category>
typename vector_traits<typename geometry_traits<Geometry>::vector_type>::direction gk_impl_direction_of(
		const Geometry& x, Category);

template<typename Geometry>
direction<geometry_traits<Geometry>::Dimension> direction_of(
		const Geometry& x) {
	return gk_impl_direction_of(x,
			geometry_traits<Geometry>::geometry_category());
}

template<typename Line>
direction<geometry_traits<Line>::Dimension> gk_impl_direction_of(const Line& l,
		line_tag) {
	typedef typename vector_traits<typename geometry_traits<Line>::vector_type>::value_type value_type;
	return normalize(l(value_type(GK_FLOAT_ONE)) - l(value_type(GK_FLOAT_ZERO)));
}

template<typename Plane>
direction<geometry_traits<Plane>::Dimension> gk_impl_direction_of(
		const Plane& plane, plane_tag) {
	typedef typename geometry_traits<Plane>::vector_type vector_type;
	typedef typename vector_traits<vector_type>::value_type value_type;
	typedef typename vector_traits<vector_type>::direction direction;

	const value_type zero = value_type(GK_FLOAT_ZERO);
	const value_type unit = value_type(GK_FLOAT_ONE);
	const direction u = normalize(plane(unit, zero));
	const direction v = normalize(plane(zero, unit));
	return cross<direction, direction, direction>()(u, v);
}

} // namespace gktypename vector_traits<typename geometry_traits<Line>::vector_type>::((

#endif /* INCLUDE_ALGORITHM_DIRECTIONOF_H_ */
