/*
 * algebra.h
 *
 *  Created on: 2016/07/29
 *      Author: J0115775
 */

#ifndef INCLUDE_VECTOR_ALGEBRA_H_
#define INCLUDE_VECTOR_ALGEBRA_H_

namespace gk {

template<std::size_t Dimension, typename T>
class vector {
public:
	vector();

	vector(const vector& other);

	template<typename Element>
	vector(const Element& element);

	~vector();

	const T& operator[](std::size_t n) const;
	T& operator[](std::size_t n);

	vector& operator=(const vector& u);

};

template<std::size_t Dimension, typename T>
vector<Dimension, T> operator+(const vector<Dimension, T>& u,
		const vector<Dimension, T>& v);
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
	matrix& operator+=(const matrix& A);
	matrix& operator-=(const matrix& A);

	template<typename Scalar>
	matrix& operator*=(const Scalar& alpha);

	template<typename Scalar>
	matrix& operator/=(const Scalar& alpha);
};

}  // namespace gk

#endif /* INCLUDE_VECTOR_ALGEBRA_H_ */
