/*
 * algorithm.h
 *
 *  Created on: 2015/10/14
 *      Author: makitaku
 */

#ifndef BSPLINE_ALGORITHM_H_
#define BSPLINE_ALGORITHM_H_

#include "bspline.h"

namespace gk {

template<typename Vector>
struct intersect_result<bspline<Vector>, bspline<Vector> > {
	typedef Vector value_type;
};

/**
 * @brief Computes intersection points of 2 B-splines.
 *
 * @param a A B-spline.
 * @param b The other B-spline.
 * @param epsilon Tolerance.
 * @param result
 * @return Next position of last intersect element.
 */
template<typename Vector, typename Tolerance, typename OutputIterator>
OutputIterator intersect(const bspline<Vector>& a, const bspline<Vector>& b,
		const Tolerance& epsilon, OutputIterator result) {
	typedef curve_traits<bspline<Vector> > traits;
	const typename traits::boundary_type bound_a = boundary(a);
	const typename traits::boundary_type bound_b = boundary(b);

	if (!is_intersect(bound_a, bound_b, epsilon)) {
		return result;
	}

	bspline<Vector> src_a = a;
	bspline<Vector> src_b = b;

	typedef typename curve_traits<bspline<Vector> >::distance_type distance_type;
	distance_type max_a;
	distance_type max_b;
	const segment<Vector> x = linear(a, max_a);
	const segment<Vector> y = linear(b, max_b);

	if (max_a < epsilon && max_b < epsilon) {
		result = intersect(x, y, epsilon, result);

	} else {
		if (max_a < epsilon) {
			result = intersect(x, src_b.subdivide(GK::Lower), epsilon, result);
			result = intersect(x, src_b, epsilon, result);

		} else if (max_b < epsilon) {
			result = intersect(src_a.subdivide(GK::Lower), y, epsilon, result);
			result = intersect(src_a, y, epsilon, result);

		} else {
			const bspline<Vector> upper_a = src_a.subdivide(GK::Lower);
			const bspline<Vector> upper_b = src_b.subdivide(GK::Lower);

			result = intersect(upper_a, upper_b, epsilon, result);
			result = intersect(upper_a, src_b, epsilon, result);
			result = intersect(src_a, upper_b, epsilon, result);
			result = intersect(src_a, src_b, epsilon, result);
		}
	}

	return result;
}

template<typename Vector, typename Tolerance, typename OutputIterator>
OutputIterator intersect(const bspline<Vector>& a, const line<Vector>& b,
		const Tolerance& epsilon, OutputIterator result) {
	return result;
}

template<typename Vector, typename Tolerance, typename OutputIterator>
OutputIterator intesect(const line<Vector>& a, const bspline<Vector>& b,
		const Tolerance& epsilon, OutputIterator result) {
	return intersect(b, a, epsilon, result);
}

template<typename Vector, typename Tolerance, typename OutputIterator>
OutputIterator intersect(const bspline<Vector>& a, const segment<Vector>& b,
		const Tolerance& epsilon, OutputIterator result) {
	const typename curve_traits<bspline<Vector> >::boundary_type bound_a =
			boundary(a);
	const typename curve_traits<segment<Vector> >::boundary_type bound_b =
			boundary(b);

	if (!is_intersect(bound_a, bound_b, epsilon)) {
		return result;
	}

	typename curve_traits<bspline<Vector> >::distance_type max_distance;
	const segment<Vector> x = linear(a, max_distance);
	if (max_distance < epsilon) {
		result = intersect(x, b, epsilon, result);
	} else {
		bspline<Vector> src_a = a;
		result = intersect(src_a.subdivide(GK::Lower), b, epsilon, result);
		result = intersect(src_a, b, epsilon, result);
	}

	return result;
}

template<typename Vector, typename Tolerance, typename OutputIterator>
OutputIterator intersect(const segment<Vector>& a, const bspline<Vector>& b,
		const Tolerance& epsilon, OutputIterator result) {
	return intersect(b, a, epsilon, result);
}

}  // namespace gk

#endif /* BSPLINE_ALGORITHM_H_ */
