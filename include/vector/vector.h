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

template<typename T, std::size_t DimensionSize>
struct vector_type {
	typedef Eigen::Matrix<T, 1, DimensionSize, Eigen::RowMajor> type;

	typedef T value_type;
	static const std::size_t Dimension = DimensionSize;
};

template<typename Vector>
struct vector_traits {
	typedef typename Vector::value_type value_type;
	static const std::size_t Dimension = Vector::Dimension;
};

template<typename T, std::size_t DimensionSize>
struct vector_traits<typename vector_type<T, DimensionSize>::type> {
	typedef T value_type;
	static const std::size_t Dimension = DimensionSize;
};

template<typename T, std::size_t Dimension>
T norm(const typename vector_type<T, Dimension>::type& v);

template<typename S, typename T, std::size_t Dimension>
typename multiplies_result<S, T>::value_type dot(
		const typename vector_type<S, Dimension>::type& u,
		const typename vector_type<T, Dimension>::type& v);

template<std::size_t Dimension>
class direction {
public:
	typedef float_type value_type;
	typedef typename vector_type<value_type, Dimension>::type vector_type;

private:
	template<typename T>
	static vector_type Normalized_(const vector_type& x) {

	}

	template<typename Vector>
	static vector_type Normalized_(const Vector& v) {
		typedef typename vector_traits<Vector>::value_type value_type;
		const typename multiplies_result<value_type, value_type>::value_type L2 =
				dot(v, v);

		const value_type L = std::sqrt(L2);

		return v / L;
	}

public:
	direction() :
			x_() {
	}

	direction(const direction& other) :
			x_(other.x_) {
	}

	template<typename T>
	direction(const typename vector_type<T, Dimension>::type& v) :
			x_() {
		const typename multiplies_result<T, T>::value_type L2 = this->x_.dot(
				this->x_);

		const T L = std::sqrt(L2);

		this->x_ = this->x_ / L;
	}

	~direction() {
	}

private:
	vector_type x_;
};

template<std::size_t RowSize, std::size_t ColumnSize, typename T>
class matrix {
public:
	matrix();
	matrix(const matrix& other);
	~matrix();

	vector_type operator()(std::size_t row) const;
	void operator()(std::size_t row, const vector_type v);

	const T& operator()(std::size_t row, std::size_t column) const;
	T& operator()(std::size_t row, std::size_t column);

	matrix& operator=(const matrix& A);

	template<typename Element>
	matrix& operator=(const Element& element);

	matrix& operator+=(const matrix& A);
	matrix& operator-=(const matrix& A);

	template<typename Scalar>
	matrix& operator*=(const Scalar& alpha);

	template<typename Scalar>
	matrix& operator/=(const Scalar& alpha);
};

}  // namespace gk

/*
 * Implements
 */

#endif /* INCLUDE_VECTOR_ALGEBRA_H_ */
