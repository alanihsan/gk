/*
 * algorithm.h
 *
 *  Created on: 2015/10/14
 *      Author: makitaku
 */

#ifndef BSPLINE_ALGORITHM_H_
#define BSPLINE_ALGORITHM_H_

#include "bspline.h"
#include "../primitive/line.h"
#include "../algorithm/kernel.h"
#include "../gkintersect.h"

namespace gk {

template<typename Vector, typename Parameter>
Parameter nearest(const bspline<Vector, Parameter>& r, const Vector& v) {
	const aabb<Vector> box = boundary(r);

}

template<typename Vector, typename Parameter, typename OutputIterator>
OutputIterator position_at(const bspline<Vector, Parameter>& r,
		OutputIterator result) {

}

/*
 * Intersection Algorithm for B-spline.
 */

template<typename Vector, typename KnotVector1, typename KnotVector2>
struct intersect_result<bspline<Vector, KnotVector1>,
		bspline<Vector, KnotVector2> > {
	typedef Vector value_type;
};

namespace impl {

/**
 * @brief Computes intersections of a B-spline and a segment.
 * @param a
 * @param box_pt1
 * @param box_pt2
 * @param epsilon
 * @param result
 * @return
 */
template<typename Vector, typename Parameter, typename Tolerance,
		typename OutputIterator>
OutputIterator intersect_bspline_segment(const bspline<Vector, Parameter>& a,
		const Vector& box_pt1, const Vector& box_pt2, const Tolerance& epsilon,
		OutputIterator result) {
	const typename curve_traits<bspline<Vector, Parameter> >::boundary_type bound_a =
			boundary(a);
	const aabb<Vector> segment_box = make_boundary(box_pt1, box_pt2);

	if (!is_intersect(bound_a, segment_box, epsilon)) {
		return result;
	}

	typename vector_traits<Vector>::value_type max_distance;
	const std::pair<Vector, Vector> X = linearize(a, max_distance);
	if (max_distance < epsilon) {
		*result = alg::intersect_2segments(X.first, X.second, box_pt1, box_pt2,
				epsilon, result);
	} else {
		bspline<Vector, Parameter> src_a = a;
		result = intersect_bspline_segment(src_a.subdivide(GK::Lower), box_pt1,
				box_pt2, epsilon, result);
		result = intersect_bspline_segment(src_a, box_pt1, box_pt2, epsilon,
				result);
	}

	return result;
}

template<typename Vector, typename KnotVector, typename Line,
		typename Tolerance, typename OutputIterator>
OutputIterator intersect_kernel(const bspline<Vector, KnotVector>& a,
		const Line& b, const Tolerance& epsilon, OutputIterator result,
		line_tag) {
	typedef typename vector_traits<Vector>::value_type value_type;

	const aabb<Vector> a_box = boundary(a);
	const Vector ref = b(value_type(GK_FLOAT_ZERO));
	const direction<vector_traits<Vector>::Dimension> u = direction_of(b);
	std::vector<Vector> X;
	alg::intersect_line_box(ref, u, a_box.min(), a_box.max(), epsilon,
			std::inserter(X, X.begin()));

	const size_t OneIntersection = 1;
	const size_t TwoIntersections = 2;

	switch (X.size()) {
	case OneIntersection:
		return intersect_bspline_segment(a, X.front(), X.front(), epsilon,
				result);

	case TwoIntersections:
		return intersect_bspline_segment(a, X[0], X[1], epsilon, result);

	default:
		return result;
	}
}

}  // namespace impl

/**
 * @brief Computes intersection points of 2 B-splines.
 *
 * @param a A B-spline.
 * @param b The other B-spline.
 * @param epsilon Tolerance.
 * @param result
 * @return Next position of last intersect element.
 */
template<typename Vector, typename Parameter1, typename Parameter2,
		typename Tolerance, typename OutputIterator>
OutputIterator intersect(const bspline<Vector, Parameter1>& a,
		const bspline<Vector, Parameter2>& b, const Tolerance& epsilon,
		OutputIterator result) {
	const aabb<Vector> bound_a = boundary(a);
	const aabb<Vector> bound_b = boundary(b);

	if (!is_intersect(bound_a, bound_b, epsilon)) {
		return result;
	}

	bspline<Vector, Parameter1> src_a = a;
	bspline<Vector, Parameter2> src_b = b;

	typedef typename vector_traits<Vector>::value_type distance_type;
	distance_type max_a;
	distance_type max_b;

	const std::pair<Vector, Vector> X = linear(a, max_a);
	const segment<Vector> x_segment(X.first, X.second);

	const std::pair<Vector, Vector> Y = linear(b, max_b);
	const segment<Vector> y_segment(Y.first, Y.second);

	if (max_a < epsilon && max_b < epsilon) {
		result = intersect(x_segment, y_segment, epsilon, result);

	} else {
		if (max_a < epsilon) {
			result = intersect(x_segment, src_b.subdivide(GK::Lower), epsilon,
					result);
			result = intersect(x_segment, src_b, epsilon, result);

		} else if (max_b < epsilon) {
			result = intersect(src_a.subdivide(GK::Lower), y_segment, epsilon,
					result);
			result = intersect(src_a, y_segment, epsilon, result);

		} else {
			const bspline<Vector, Parameter1> upper_a = src_a.subdivide(
					GK::Lower);
			const bspline<Vector, Parameter2> upper_b = src_b.subdivide(
					GK::Lower);

			result = intersect(upper_a, upper_b, epsilon, result);
			result = intersect(upper_a, src_b, epsilon, result);
			result = intersect(src_a, upper_b, epsilon, result);
			result = intersect(src_a, src_b, epsilon, result);
		}
	}

	return result;
}

template<typename Vector, typename Parameter, typename Other,
		typename Tolerance, typename OutputIterator>
OutputIterator intersect(const bspline<Vector, Parameter>& a, const Other& b,
		const Tolerance& epsilon, OutputIterator result) {
	return impl::intersect_kernel(a, b, epsilon, result,
			typename geometry_traits<Other>::geometry_category());
}

template<typename Vector, typename Parameter, typename Other,
		typename Tolerance, typename OutputIterator>
OutputIterator intersect(const Other& a, const bspline<Vector, Parameter>& b,
		const Tolerance& epsilon, OutputIterator result) {
	OutputIterator end = intersect(b, a, epsilon, result);
	std::reverse(result, end);
	return result;
}

}  // namespace gk

#endif /* BSPLINE_ALGORITHM_H_ */
