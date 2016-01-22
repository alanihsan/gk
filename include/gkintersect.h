/*
 * gkintersect.h
 *
 *  Created on: 2015/10/06
 *      Author: makitaku
 */

#ifndef GKINTERSECT_H_
#define GKINTERSECT_H_

namespace gk {

/**
 * @brief Result type of intersection.
 *
 * @tparam Geometry1
 * @tparam Geometry2
 */
template<typename Geometry1, typename Geometry2>
struct intersect_result {
	typedef void type;
};

/**
 * @brief Template function to intersect 2 geometries.
 *
 * @param a First geometry.
 * @param b Second geometry.
 * @param epsilon Tolerance value.
 * @param result
 * @return
 */
template<typename Geometry1, typename Geometry2, typename Tolerance,
		typename OutputIterator>
OutputIterator intersect(const Geometry1& a, const Geometry2& b,
		const Tolerance& epsilon, OutputIterator result);

}  // namespace gk

#endif /* GKINTERSECT_H_ */
