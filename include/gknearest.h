/*
 * gknearest.h
 *
 *  Created on: 2016/01/22
 *      Author: makitaku
 */

#ifndef GKNEAREST_H_
#define GKNEAREST_H_

#include "gkdirection.h"

namespace gk {

namespace inner {

template<typename Geometry, typename Vector, typename GeometryCategory>
Vector nearest(const Geometry& g, const Vector& v, GeometryCategory);

template<typename Vector>
Vector nearest_to_line(const Vector& reference,
		const direction<vector_traits<Vector>::Dimension>& u, const Vector& v) {
	const Vector r = v - reference;
	const dot<Vector, direction<vector_traits<Vector>::Dimension>,
			typename vector_traits<Vector>::value_type> dot;
	return dot(r, u) * u - r;
}

template<typename Line, typename Vector>
Vector nearest(const Line& line, const Vector& v, line_tag) {
	typedef typename vector_traits<Vector>::value_type value_type;
	return nearest_to_line(line(value_type(GK_FLOAT_ZERO)), direction_of(line),
			v);
}

template<typename Segment, typename Vector>
Vector nearest(const Segment& segment, const Vector& v, segment_tag) {
	typedef typename vector_traits<Vector>::value_type value_type;
	typedef typename direction<vector_traits<Vector>::Dimension> direction;

	const direction u = direction_of(segment);
	const Vector r = nearest_to_line(segment[GK::StartEdge], u, v);

	const dot<Vector, direction, value_type> dot;
	const value_type t = dot(r, u);
	if (std::signbit(t)) {
		return segment[GK::StartEdge] - v;
	} else {
		return (t >= length(segment)) ? segment[GK::EndEdge] - v : r;
	}
}

template<typename Plane, typename Vector>
Vector nearest(const Plane& plane, const Vector& v, plane_tag) {
	typedef typename vector_traits<Vector>::value_type value_type;
	typedef direction<vector_traits<Vector>::Dimension> direction;

	const direction n = direction_of(plane);
	const Vector r = plane(value_type(GK_FLOAT_ZERO), value_type(GK_FLOAT_ZERO))
			- v;

	const dot<Vector, direction, value_type> dot;
	return dot(r, n) * n;
}

}  // namespace inner

/**
 * @brief Compute a vector being from a position vector @a v to a nearest position
 * on a geometry @a g.
 * @param g
 * @param v
 * @return
 */
template<typename Geometry, typename Vector>
Vector nearest(const Geometry& g, const Vector& v) {
	return inner::nearest(g, v, geometry_traits<Geometry>::geometry_category());
}

namespace inner {

template<typename Geometry1, typename Geometry2, typename Geometry1Category,
		typename Geometry2Category>
std::pair<
		typename vector_traits<typename geometry_traits<Geometry1>::vector_type>::value_type,
		typename vector_traits<typename geometry_traits<Geometry2>::vector_type>::value_type> nearest_between(
		const Geometry1& a, const Geometry2& b, Geometry1Category,
		Geometry2Category);

template<typename Vector>
std::pair<typename vector_traits<Vector>::value_type,
		typename vector_traits<Vector>::value_type> nearest_between_lines(
		const Vector& reference1,
		const direction<vector_traits<Vector>::Dimension>& direction1,
		const Vector& reference2,
		const direction<vector_traits<Vector>::Dimension>& direction2) {

	typedef typename vector_traits<Vector>::value_type length;
	typedef direction<vector_traits<Vector>::Dimension> direction;

	const gkfloat alpha = dot<direction, direction, gkfloat>(direction1,
			direction2);
	const gkfloat beta = GK_FLOAT_ONE - alpha * alpha;

	const Vector r = reference1 - reference2;

	if (beta == GK_FLOAT_ZERO) {
		/* parallel */
		return std::make_pair(length(GK_FLOAT_ZERO), dot(r, direction2));

	} else {
		/* no parallel */
		const dot<Vector, direction, length> dot;
		const length s = (alpha * dot(r, direction2) - dot(r, direction1))
				/ beta;
		const length t = dot(r, direction2) + alpha * s;

		return std::make_pair(s, t);
	}
}

template<typename Line>
std::pair<
		typename vector_traits<typename geometry_traits<Line>::vector_type>::value_type,
		typename vector_traits<typename geometry_traits<Line>::vector_type>::value_type> nearest_between(
		const Line& l, const Line& m, line_tag, line_tag) {
	typedef typename geometry_traits<Line>::vector_type vector_type;
	typedef typename vector_traits<vector_type>::value_type length;
	typedef direction<vector_traits<vector_type>::Dimension> direction;

	return nearest_between_lines(l(length(GK_FLOAT_ZERO)), direction_of(l),
			m(length(GK_FLOAT_ZERO)), direction_of(m));
}

template<typename Segment>
std::pair<typename geometry_traits<Segment>::vector_type,
		typename geometry_traits<Segment>::vector_type> nearest_between(
		const Segment& l, const Segment& m, segment_tag, segment_tag) {
	typedef typename geometry_traits<Segment>::vector_type vector_type;
	typedef typename vector_traits<vector_type>::value_type length;

	const std::pair<length, length> r = nearest_between_lines(
			l(length(GK_FLOAT_ZERO)), direction_of(l), m(length(GK_FLOAT_ZERO)),
			direction_of(m));

	const length s = (std::signbit(r.first)) ? length(GK_FLOAT_ZERO) :
						(r.first >= length(l)) ? length(l) : r.first;

	const length t = (std::signbit(r.second)) ? length(GK_FLOAT_ZERO) :
						(r.second >= length(m)) ? length(m) : r.second;

	return std::make_pair(s, t);
}

}  // namespace inner

/**
 * @brief Compute
 * @param a
 * @param b
 * @return
 */
template<typename Geometry1, typename Geometry2>
std::pair<typename geometry_traits<Geometry1>::vector_type,
		typename geometry_traits<Geometry2>::vector_type> nearest_betweeen(
		const Geometry1& a, const Geometry2& b) {
	return inner::nearest_between(a, b,
			geometry_traits<Geometry1>::geometry_category(),
			geometry_traits<Geometry2>::geometry_category());
}

}  // namespace gk

#endif /* GKNEAREST_H_ */
