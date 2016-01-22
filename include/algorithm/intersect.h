/*
 * intersect.h
 *
 *  Created on: 2015/07/15
 *      Author: tmakimoto
 */

#ifndef ALGORITHM_INTERSECT_H_
#define ALGORITHM_INTERSECT_H_

//#include "../gkscalar.h"
#include "../gkdirection.h"
#include "../gkaabb.h"

#include <utility>

namespace gk {

/**
 * @brief Computes the intersection of 2 lines.
 *
 * Each line is represented from its edge position vectors.
 *
 * @param start1 Start position vector of the first segment.
 * @param end1 End position vector of the first segment.
 * @param start2 Start position vector of the second segment.
 * @param end2 End position vector of the second segment.
 * @param epsilon
 * @return
 */
template<typename Vector, typename Tolerance>
Vector intersect(const Vector& start1, const Vector& end1, const Vector& start2,
		const Vector& end2,
		const Tolerance& epsilon = std::numeric_limits<
				typename vector_traits<Vector>::value_type>::epsilon()) { const
	Vector u = end1 - start1;
	const Vector v = end2 - start2;
}

template<typename Vector>
std::pair<bool, Vector> intersect(const std::pair<Vector, Vector>& l,
		const std::pair<Vector, Vector>& m);

namespace impl {

template<typename Line>
bool is_intersect_line_aabb_impl(const Line& line,
		const aabb<typename geometry_traits<Line>::vector_type>& box,
		dimension<GK::GK_2D>) {
	typedef typename geometry_traits<Line>::vector_type vector_type;
	typedef typename scalar_traits<
			typename vector_traits<vector_type>::value_type>::float_type float_type;

	const vector_type r = line(float_type(GK_FLOAT_ZERO));
	const vector_type v[] = { box.min() - r, vector_type(box.max()[GK::X],
			box.min()[GK::Y]) - r, box.max() - r, vector_type(box.min()[GK::X],
			box.max()[GK::Y]) - r };

	const direction<vector_type> d = get_direction(line);
	typename cross_result<vector_type, direction<vector_type> >::type n[] = {
			cross(v[0], d), cross(v[1], d), cross(v[2], d), cross(v[3], d) };
	return std::signbit(n[0]) ^ std::signbit(n[1]) ^ std::signbit(n[2])
			^ std::signbit(n[3]);
}

template<typename Line>
bool is_intersect_line_aabb_impl(const Line& line,
		const aabb<typename geometry_traits<Line>::vector_type>& box,
		dimension<GK::GK_3D>) {

	typedef typename geometry_traits<Line>::vector_type vector_type;
	typedef typename vector_traits<vector_type>::value_type value_type;
	typedef typename scalar_traits<value_type>::float_type float_type;

	const direction<vector_type> t = get_direction(line);

	const vector_type r = line(vector_type(GK_FLOAT_ZERO));
	const vector_type min_vector = box.min();
	const vector_type max_vector = box.max();

	value_type input[GK::GK_3D] = { value_type(GK_FLOAT_ZERO) };
	value_type output[GK::GK_3D] =
			{
					value_type(
							std::numeric_limits<
									typename scalar_traits<value_type>::float_type>::infinity()) };

	{
		const value_type x_min = (r[GK::X] - min_vector[GK::X]) / t[GK::X];
		const value_type x_max = (r[GK::X] - max_vector[GK::X]) / t[GK::X];
		(x_min < x_max) ?
				input[GK::X] = x_min, output[GK::X] = x_max : input[GK::X] =
						x_max, output[GK::X] = x_min;
	}

	{
		const value_type y_min = (r[GK::Y] - min_vector[GK::Y]) / t[GK::Y];
		const value_type y_max = (r[GK::Y] - max_vector[GK::Y]) / t[GK::Y];
		(y_min < y_max) ?
				input[GK::Y] = y_min, output[GK::Y] = y_max : input[GK::Y] =
						y_max, output[GK::Y] = y_min;
	}

	{
		const value_type z_min = (r[GK::Z] - min_vector[GK::Z]) / t[GK::Z];
		const value_type z_max = (r[GK::Z] - max_vector[GK::Z]) / t[GK::Z];
		(z_min < z_max) ?
				input[GK::Z] = z_min, output[GK::Z] = z_max : input[GK::Z] =
						z_max, output[GK::Z] = z_min;
	}

	const value_type* max_input = std::max_element(input, input + GK::GK_3D);
	const value_type* min_output = std::min_element(output, output + GK::GK_3D);
	return (*min_output - *max_input) >= value_type(GK_FLOAT_ZERO);
}

template<typename Plane>
bool is_intersect_plane_aabb_impl(const Plane& plane,
		const aabb<typename geometry_traits<Plane>::vector_type>& aabb) {
}

template<typename Geometry1, typename Geometry2, typename Geometry1Category,
		typename Geometry2Category>
struct intersect_result_impl;

template<typename Geometry1, typename Geometry2, typename Geometry1Category,
		typename Geometry2Category, typename Dimension>
typename intersect_result_impl<Geometry1, Geometry2, Geometry1Category,
		Geometry2Category> intersect_impl(const Geometry1& a,
		const Geometry2& b, Geometry1Category, Geometry2Category, Dimension);

template<typename Line>
struct intersect_result_impl<Line, Line, line_tag, line_tag> {
	typedef typename geometry_traits<Line>::vector_type type;
};

#define GK_INTERSECT_RESULT_LINES_FORM(line1, line1_category, line2, line2_category) \
	template<typename line1, typename line2> \
	struct intersect_result_impl<line1, line2, line1_category, line2_category> { \
		typedef typename requirement<typename geometry_traits<line1>::vector_type, \
				typename geometry_traits<line2>::vector_type>::type type; \
	};

template<typename Line>
struct intersect_result_impl<Line, Line, line_tag, segment_tag> {
	typedef typename geometry_traits<Line>::vector_type type;
};

template<typename Plane>
struct intersect_result_impl<Plane, Plane, plane_tag, plane_tag> {
// The result of intersection two planes is a line.
//	typedef intersect_result_plane_plane_impl<
//			typename geometry_traits<Plane>::vector_type> value_type;
	typedef typename geometry_traits<Plane>::vector_type vector_type;
	typedef std::pair<vector_type, vector_type> type;
};

template<typename Line, typename Plane>
struct intersect_result_impl<Line, Plane, line_tag, plane_tag> {
	typedef typename requirement<typename geometry_traits<Line>::vector_type,
			typename geometry_traits<Plane>::vector_type>::type type;
};

template<typename Line>
typename intersect_result_impl<Line, Line, line_tag, line_tag>::type intersect_impl(
		const Line& a, const Line& b, line_tag, line_tag,
		dimension<GK::GK_2D>) {
	typedef typename geometry_traits<Line>::vector_type vector_type;

	const direction<vector_type> u = get_direction(a);
	const direction<vector_type> v = get_direction(b);

	const typename vector_traits<vector_type>::value_type One(GK_FLOAT_ONE);
	const direction<vector_type> w(
			vector_type(-One * v[GK::Y, One * v[GK::X]]));

	typedef typename curve_traits<Line>::parameter parameter;
	const typename curve_traits<Line>::parameter Zero(GK_FLOAT_ZERO);
	const parameter s = dot(b(Zero) - a(Zero), w) / dot(u, w);

	return a(s);
}

template<typename Line>
typename geometry_traits<Line>::vector_type intersect_impl(const Line& a,
		const Line& b, line_tag, line_tag, dimension<GK::GK_3D>) {
	typedef typename geometry_traits<Line>::vector_type vector_type;

	const direction<vector_type> u = get_direction(a);
	const direction<vector_type> v = get_direction(b);

	const direction<vector_type> w = normalize(
			vector_traits<vector_type>::value_type(
			GK_FLOAT_ONE) * cross(u, v));

	typedef typename curve_traits<Line>::parameter parameter;
	const parameter Zero(GK_FLOAT_ZERO);

	const vector_type r = a(Zero) - b(Zero);
	const vector_type p = dot(r, w) * w - r;

	const vector_type f = a(dot(p, u));
	const vector_type g = b(dot(p, v));

	const parameter Infinity(
			std::numeric_limits<typename scalar_traits<parameter>::float_type>::infinity());

	return (f == g) ? f : vector_type(Infinity, Infinity, Infinity);
}

template<typename Plane>
typename intersect_result_impl<Plane, Plane, plane_tag, plane_tag>::type intersect_impl(
		const Plane& a, const Plane& b, plane_tag, plane_tag,
		dimension<GK::GK_3D>) {
	typedef typename geometry_traits<Plane>::vector_type vector_type;

	const direction<vector_type> m = direction_of(a);
	const direction<vector_type> n = direction_of(b);

	if (m == n) {
		return std::make_pair(vector_type(), vector_type());
	}

	const vector_type u = vector_traits<vector_type>::value_type(GK_FLOAT_ONE)
			* cross(m, n);

	const typename surface_traits<Plane>::parameter zero(GK_FLOAT_ZERO);

	const vector_type r = b(zero, zero)
			+ dot(a(zero, zero) - b(zero, zero), m) * m;

	return std::make_pair(r, r + u);
}

template<typename Vector, typename Tolerance>
std::pair<Vector, bool> intersect_segment_segment_impl(const Vector& a1,
		const Vector& a2, const Vector& b1, const Vector& b2,
		const Tolerance& epsilon, dimension<GK::GK_2D>) {
	typedef typename vector_traits<Vector>::value_type value_type;
	typedef division_result<value_type, value_type> float_type;
//	const direction<Vector> u = direction_of(a2 - a1);
//	const direction<Vector> v = direction_of(b2 - b1);
//
//	const typename vector_traits<Vector>::value_type InfinityValue(
//			std::numeric_limits<typename vector_traits<Vector>::value_type>::infinity());
//	if (u == v) {
//		return std::make_pair(Vector(), false);
//	}

	std::pair<Vector, Vector> nearests = nearest_lines(a1, a2, b1, b2);
	if (norm(nearests.first - nearests.second) < epsilon) {
		const float_type s = norm(nearests.first - a1) / (norm(a2 - a1));
		const float_type t = norm(nearests.second - b1) / (norm(b2 - b1));

	} else {
		return std::make_pair(Vector(), false);
	}
}

template<typename Vector, typename Tolerance>
std::pair<Vector, bool> intersect_segment_segment_impl(const Vector& a1,
		const Vector& a2, const Vector& b1, const Vector& b2,
		const Tolerance& epsilon, dimension<GK::GK_3D>) {

}

} // namespace impl

template<typename Geometry1, typename Geometry2>
struct intersect_result {
	typedef typename impl::intersect_result_impl<Geometry1, Geometry2,
			typename geometry_traits<Geometry1>::geometry_category,
			typename geometry_traits<Geometry2>::geometry_category>::type type;
};

/**
 * @brief Computes
 * @param a
 * @param b
 * @return
 */
template<typename Geometry1, typename Geometry2>
typename intersect_result<Geometry1, Geometry2>::type intersect(
		const Geometry1& a, const Geometry2& b) {
	return impl::intersect_impl(a, b,
			geometry_traits<Geometry1>::geometry_category(),
			geometry_traits<Geometry2>::geometry_category(),
			requirement<
					typename vector_traits<
							typename geometry_traits<Geometry1>::vector_type>::dimension_category,
					typename vector_traits<
							typename geometry_traits<Geometry2>::vector_type>::dimension_category>::value_type());
}

template<typename Vector, typename Tolerance>
std::pair<Vector, bool> intersect_segment_segment(const Vector& a1,
		const Vector& a2, const Vector& b1, const Vector& b2,
		const Tolerance& epsilon = std::numeric_limits<
				vector_traits<Vector>::value_type>::epsilon()) {
	return impl::intersect_segment_segment_impl(a1, a2, b1, b2, epsilon,
			dimension<vector_traits<Vector>::Dimension>());
}

}  // namespace gk

#endif /* ALGORITHM_INTERSECT_H_ */
