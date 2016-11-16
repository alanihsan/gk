/*
 * algebra.h
 *
 *  Created on: 2016/07/29
 *      Author: J0115775
 */

#ifndef INCLUDE_VECTOR_ALGEBRA_H_
#define INCLUDE_VECTOR_ALGEBRA_H_

#include <Eigen/Core>

namespace gk {

<<<<<<< HEAD
//template<typename T, std::size_t DimensionSize>
//struct vector_type {
//	typedef Eigen::Matrix<T, 1, DimensionSize, Eigen::RowMajor> type;
//
//	typedef T value_type;
//	static const std::size_t Dimension = DimensionSize;
//};
=======
template<typename T, std::size_t DimensionSize>
struct vector_type {
	typedef Eigen::Matrix<T, 1, DimensionSize, Eigen::RowMajor> type;

	typedef T value_type;
	static const std::size_t Dimension = DimensionSize;
};
>>>>>>> 5ed5cf2db97b02bd867d82a1513c790b646cb925

template<typename Vector>
struct vector_traits {
	typedef typename Vector::value_type value_type;
	static const std::size_t Dimension = Vector::Dimension;
};

template<typename T, std::size_t DimensionSize>
<<<<<<< HEAD
struct vector_traits<Eigen::Matrix<T, DimensionSize, Eigen::RowMajor> > {
=======
struct vector_traits<typename vector_type<T, DimensionSize>::type> {
>>>>>>> 5ed5cf2db97b02bd867d82a1513c790b646cb925
	typedef T value_type;
	static const std::size_t Dimension = DimensionSize;
};

template<typename Vector, typename T, std::size_t Dimension>
void assign(const Vector& src,
<<<<<<< HEAD
		const Eigen::Matrix<T, Dimension, Eigen::RowMajor>& dst) {
=======
		const typename vector_type<T, Dimension>::type& dst) {
>>>>>>> 5ed5cf2db97b02bd867d82a1513c790b646cb925
	assign_dimension(src, dst, dimension_tag<Dimension>());
}

template<typename Vector, typename T>
void assign_dimension(const Vector& src,
<<<<<<< HEAD
		const Eigen::Matrix<T, GK::GK_2D, Eigen::RowMajor>& dst,
=======
		const typename vector_type<T, GK::GK_2D>::type& dst,
>>>>>>> 5ed5cf2db97b02bd867d82a1513c790b646cb925
		dimension_tag<GK::GK_2D>) {
	dst[GK::X] = src[GK::X];
	dst[GK::Y] = src[GK::Y];
}

template<typename Vector, typename T>
void assign_dimension(const Vector& src,
<<<<<<< HEAD
		const Eigen::Matrix<T, GK::GK_3D, Eigen::RowMajor>& dst,
=======
		const typename vector_type<T, GK::GK_3D>::type& dst,
>>>>>>> 5ed5cf2db97b02bd867d82a1513c790b646cb925
		dimension_tag<GK::GK_3D>) {
	dst[GK::X] = src[GK::X];
	dst[GK::Y] = src[GK::Y];
	dst[GK::Z] = src[GK::Z];
}

template<typename T, std::size_t Dimension>
<<<<<<< HEAD
T norm(const Eigen::RowMajor<T, Dimension, Eigen::RowMajor>& v);

template<typename S, typename T, std::size_t Dimension>
typename multiplies_result<S, T>::value_type dot(
		const Eigen::Matrix<S, Dimension, Eigen::RowMajor>& u,
		const Eigen::Matrix<T, Dimension, Eigen::RowMajor>& v);
=======
T norm(const typename vector_type<T, Dimension>::type& v);

template<typename S, typename T, std::size_t Dimension>
typename multiplies_result<S, T>::value_type dot(
		const typename vector_type<S, Dimension>::type& u,
		const typename vector_type<T, Dimension>::type& v);
>>>>>>> 5ed5cf2db97b02bd867d82a1513c790b646cb925

template<std::size_t Dimension>
class direction {
public:
	typedef float_type value_type;
<<<<<<< HEAD
	typedef Eigen::Matrix<value_type, Dimension, Eigen::RowMajor> vector_type;
=======
	typedef typename vector_type<value_type, Dimension>::type vector_type;
>>>>>>> 5ed5cf2db97b02bd867d82a1513c790b646cb925

private:
	template<typename T>
	static vector_type Normalized_(const vector_type& x) {

	}

	template<typename Vector>
	static vector_type Normalized_(const Vector& v) {
		typedef typename vector_traits<Vector>::value_type T;
		const typename multiplies_result<T, T>::value_type L2 = dot(v, v);

		const T L = std::sqrt(L2);

		return v / L;
	}

public:
	direction() :
			x_() {
	}

	direction(const direction& other) :
			x_(other.x_) {
	}

	template<typename Vector>
	direction(const Vector& v) :
			x_(Normalized_(v)) {
	}

	~direction() {
	}

	value_type operator[](std::size_t n) const {
		return this->x_[n];
	}

	value_type& operator[](std::size_t n) {
		return this->x_[n];
	}

	direction& operator=(const direction& rhs) {
		if (&rhs == this) {
			return *this;
		}

		this->x_ = rhs.x_;
		return *this;
	}

private:
	vector_type x_;
};

template<typename T, std::size_t Dimension>
<<<<<<< HEAD
Eigen::Matrix<T, Dimension, Eigen::RowMajor> operator*(const T& alpha,
		const direction<Dimension>& d) {
	Eigen::Matrix<T, Dimension, Eigen::RowMajor> vector_t;
	vector_t r;
	for (std::size_t i = 0; i < Dimension; ++i) {
		r[i] = alpha * d[i];
	}
=======
typename vector_type<T, Dimension>::type operator*(const T& alpha,
		const direction<Dimension>& d) {
	typedef typename vector_type<T, Dimension>::type vector_t;
	vector_t r;
	for (std::size_t i = 0; i < Dimension; ++i)
		r[i] = alpha * d[i];
>>>>>>> 5ed5cf2db97b02bd867d82a1513c790b646cb925

	return r;
}

template<typename T, std::size_t Dimension>
<<<<<<< HEAD
Eigen::Matrix<T, Dimension, Eigen::RowMajor> operator*(
=======
typename vector_type<T, Dimension>::type operator*(
>>>>>>> 5ed5cf2db97b02bd867d82a1513c790b646cb925
		const direction<Dimension>& d, const T& alpha) {
	return alpha * d;
}

}  // namespace gk

#endif /* INCLUDE_VECTOR_ALGEBRA_H_ */
