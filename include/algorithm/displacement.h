/*
 * displacement.h
 *
 *  Created on: 2015/07/07
 *      Author: tmakimoto
 */

#ifndef ALGORITHM_DISPLACEMENT_H_
#define ALGORITHM_DISPLACEMENT_H_

//#include "../gkprimitive.h"

namespace gk {

namespace impl {

template<typename Geometry1, typename Geometry2>
struct displacement_result_impl {
	typedef typename requirement<
			typename geometry_traits<Geometry1>::vector_type,
			typename geometry_traits<Geometry2>::vector_type>::value_type type;
};

template<typename Geometry1, typename Geometry2, typename Geometry1Category,
		typename Geometry2Category>
typename displacement_result_impl<Geometry1, Geometry2>::type diplacement_impl(
		const Geometry1& first, const Geometry2& second, Geometry1Category,
		Geometry2Category);

template<typename Line>
typename displacement_result_impl<Line, Line>::type displacement_line_line_impl(
		const Line& l, const Line& m, dimension<GK::GK_2D>) {
	return curve_traits<Line>::vector_type(); // returns zero vector.
}

template<typename Line>
typename displacement_result_impl<Line, Line>::type displacement_line_line_impl(
		const Line& l, const Line& m, dimension<GK::GK_3D>) {

}

template<typename Line>
typename displacement_result_impl<Line, Line>::type diplacement_impl(
		const Line& first, const Line& second, line_tag, line_tag) {
	return displacement_line_line_impl(first, second,
			vector_traits<typename geometry_traits<Line>::vector_type>::dimension_category());
}

} // namespace impl

template<typename Geometry1, typename Geometry2>
struct displacement_result {
	typedef typename impl::displacement_result_impl<Geometry1, Geometry2>::type type;
};

/**
 * @brief Computes the displacement from @a first to @a second.
 * @param first
 * @param second
 * @return
 */
template<typename Geometry1, typename Geometry2>
typename displacement_result<Geometry1, Geometry2>::type displacement(
		const Geometry1& first, const Geometry2& second);

template<typename Vector>
Vector displacement(const Vector& position, const Vector& segment_start,
		const Vector&segment_end) {

}

} // namespace gk

#endif /* ALGORITHM_DISPLACEMENT_H_ */
