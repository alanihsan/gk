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

template<std::size_t Dimension, typename T>
class vector {
private:
	Eigen::Matrix<1, Dimension, Eigen::RowMajor> type;

public:
	vector() :
			x_() {
	}

	vector(const vector& other) :
			x_(other.x_) {
	}

	template<typename Element>
	vector(const Element& v) :
			x_() {
		this->assign_(v, dimension_tag<Dimension>());
	}

	~vector() {
	}

	const T& operator[](std::size_t n) const {
		return this->x_[n];
	}

	T& operator[](std::size_t n) {
		return this->x_[n];
	}

	vector& operator=(const vector& v) {
		if (&v == this) {
			return *this;
		}

		this->x_ = v.x_;
		return *this;
	}

	template<typename Element>
	vector& operator=(const Element& v) {
		this->assign_(v, dimension_tag<Dimension>());
	}

	vector& operator+=(const vector& v) {
		this->x_ += v.x_;
		return *this;
	}

	vector& operator-=(const vector& v) {
		this->x_ -= v.x_;
		return *this;
	}

	template<typename Scalar>
	vector& operator*=(const Scalar& alpha) {
		this->x_ *= alpha;
		return *this;
	}

	template<typename Scalar>
	vector& operator/=(const Scalar& alpha) {
		this->x_ /= alpha;
		return *this;
	}

private:
	type x_;

private:
	template<typename Element>
	void assign_(const Element& v, dimension_tag<GK::GK_2D>) {
		this->x_[GK::X] = v[GK::X];
		this->x_[GK::Y] = v[GK::Y];
	}

	template<typename Element>
	void assign_(const Element& v, dimension_tag<GK::GK_3D>) {
		this->x_[GK::X] = v[GK::X];
		this->x_[GK::X] = v[GK::X];
		this->x_[GK::X] = v[GK::X];
	}
};

template<std::size_t Dimension, typename T>
vector<Dimension, T> operator+(const vector<Dimension, T>& u,
		const vector<Dimension, T>& v) {
	vector<Dimension, T> r = u;
	r += v;
	return r;
}

template<std::size_t Dimension, typename T>
vector<Dimension, T> operator-(const vector<Dimension, T>& u,
		const vector<Dimension, T>& v);

template<std::size_t Dimension, typename T, typename Scalar>
vector<Dimension, T> operator*(const Scalar& alpha,
		const vector<Dimension, T>& v);

template<std::size_t Dimension, typename T, typename Scalar>
vector<Dimension, T> operator*(const vector<Dimension, T>& v,
		const Scalar& alpha);

template<std::size_t Dimension, typename T, typename Scalar>
vector<Dimension, T> operator/(const vector<Dimension, T>& v,
		const Scalar& alpha);

template<std::size_t Dimension, typename T>
T norm(const vector<Dimension, T>& v);

template<std::size_t Dimension, typename S, typename T>
typename multiplies_result<S, T>::value_type dot(const vector<Dimension, S>& u,
		const vector<Dimension, T>& v);

template<std::size_t RowSize, std::size_t ColumnSize, typename T>
class matrix {
public:
	matrix();
	matrix(const matrix& other);
	~matrix();

	vector operator()(std::size_t row) const;
	void operator()(std::size_t row, const vector v);

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
