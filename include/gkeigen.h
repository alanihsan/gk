/*
 * gkeigen.h
 *
 *  Created on: 2016/01/22
 *      Author: makitaku
 */

#ifndef GKEIGEN_H_
#define GKEIGEN_H_

#ifdef GK_USING_EIGEN

#include "gkvector.h"
#include <Eigen/Core>

namespace gk {

template<typename Scalar, size_t DimensionSize>
struct vector_traits<Eigen::Matrix<Scalar, 1, DimensionSize> > {
	typedef Scalar value_type;

	typedef value_type* iterator;
	typedef const value_type* const_iterator;
	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	static const size_t Dimension = DimensionSize;
	static const bool IsHomegeneous = false;
};

template<typename Scalar, size_t DimensionSize>
typename vector_traits<Eigen::Matrix<Scalar, 1, DimensionSize> >::const_iterator begin(
		const Eigen::Matrix<Scalar, 1, DimensionSize>& v) {
	return v.data();
}

template<typename Scalar, size_t DimensionSize>
typename vector_traits<Eigen::Matrix<Scalar, 1, DimensionSize> >::iterator begin(
		Eigen::Matrix<Scalar, 1, DimensionSize>& v) {
	return v.data();
}

template<typename Scalar, size_t DimensionSize>
typename vector_traits<Eigen::Matrix<Scalar, 1, DimensionSize> >::const_iterator end(
		const Eigen::Matrix<Scalar, 1, DimensionSize>& v) {
	return v.data() + DimensionSize;
}

template<typename Scalar, size_t DimensionSize>
typename vector_traits<Eigen::Matrix<Scalar, 1, DimensionSize> >::iterator end(
		Eigen::Matrix<Scalar, 1, DimensionSize>& v) {
	return v.data() + DimensionSize;
}

template<typename Scalar, size_t DimensionSize>
Scalar norm(const Eigen::Matrix<Scalar, 1, DimensionSize>& v) {
	return v.norm();
}

}  // namespace gk

#endif

#endif /* GKEIGEN_H_ */
