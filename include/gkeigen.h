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

template<size_t DimensionSize>
struct eigen_direction {
	typedef Eigen::Matrix<gkfloat, 1, DimensionSize> vector_type;
	typedef typename vector_type::Scalar value_type;
	typedef const value_type* const_pointer;
	typedef const_pointer const_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	static const size_t Dimension = DimensionSize;

	const vector_type X;

	eigen_direction(const vector_type& A) :
			X(A / A.norm()) {
	}

	~eigen_direction() {
	}

	const_iterator begin() const {
		return this->X.data();
	}

	const_iterator end() const {
		return this->X.data() + Dimension;
	}

	const_reverse_iterator rbegin() const {
		return const_reverse_iterator(this->end());
	}

	const_reverse_iterator rend() const {
		return const_reverse_iterator(this->begin());
	}

	const value_type& operator[](size_t n) const {
		return this->X[n];
	}

private:
	eigen_direction();
};

template<typename Scalar, size_t DimensionSize>
struct vector_traits<Eigen::Matrix<Scalar, 1, DimensionSize> > {
	typedef Scalar value_type;

	typedef value_type* iterator;
	typedef const value_type* const_iterator;
	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	static const size_t Dimension = DimensionSize;
	static const bool IsHomegeneous = false;

	typedef eigen_direction<DimensionSize> direction;
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

template<typename Scalar, size_t DimensionSize>
Eigen::Matrix<Scalar, 1, DimensionSize> operator*(const Scalar& alpha,
		const eigen_direction<DimensionSize>& u) {
	return alpha * u.X;
}

template<typename Scalar, size_t DimensionSize>
Eigen::Matrix<Scalar, 1, DimensionSize> operator*(
		const eigen_direction<DimensionSize>& u, const Scalar& alpha) {
	return alpha * u;
}

}  // namespace gk

#endif

#endif /* GKEIGEN_H_ */
