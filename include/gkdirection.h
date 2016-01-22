/*
 * gkdirection.h
 *
 *  Created on: 2015/10/01
 *      Author: makitaku
 */

#ifndef GKDIRECTION_H_
#define GKDIRECTION_H_

#include <gkvector.h>
#include <functional>
#include "gkfunctional.h"
#include "gkgeometry.h"

namespace gk {

/**
 * @brief
 *
 * @author Takuya Makimoto
 * @date 2015
 */
/*template<typename Vector>*/
template<size_t DimensionSize>
class direction {
public:
	typedef gkfloat value_type;
	typedef const value_type* const_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	static const size_t Dimension = DimensionSize;/*vector_traits<Vector>::Dimension;*/
	static const size_t ElementSize = DimensionSize;/*vector_traits<Vector>::Dimension; //dimension_category::Size;*/

private:

	/**
	 * @brief Normalize a vector.
	 * @param v
	 * @param d
	 */
	template<typename Vector>
	static void normalize_(const Vector& v, direction& d) {
		std::transform(begin(v), end(v), d.x_,
				std::bind2nd(
						divides<value_type,
								typename vector_traits<Vector>::value_type>(),
						value_type(GK_FLOAT_ONE) / norm(v)));
	}

public:
	/**
	 * @brief Default constructor.
	 */
	direction() :
			x_() {
		std::fill(this->x_, this->x_ + Dimension, value_type(GK_FLOAT_ZERO));
	}

	/**
	 * @brief Copy constructor.
	 * @param other
	 */
	direction(const direction& other) :
			x_() {
		std::copy(other.x_, other.x_ + ElementSize, this->x_);
	}

	template<typename Vector>
	explicit direction(const Vector& v) :
			x_() {
		normalize_(v, *this);
	}

	~direction() {
	}

	value_type operator[](const size_t n) const {
		return this->x_[n];
	}

	direction& operator=(const direction& u) {
		if (&u == this) {
			return *this;
		}

		std::copy(u.x_, u.x_ + ElementSize, this->x_);
		return *this;
	}

	const_iterator begin() const {
		return this->x_;
	}

	const_iterator end() const {
		return this->x_ + Dimension;
	}

	const_reverse_iterator rbegin() const {
		return std::reverse_iterator<const_iterator>(this->end());
	}

	const_reverse_iterator rend() const {
		return std::reverse_iterator<const_iterator>(this->begin());
	}

	const value_type* data() const {
		return this->x_;
	}

private:
//	value_type x_[ElementSize];
	gkfloat x_[ElementSize];
};

template<typename Vector>
vector<vector_traits<Vector>::Dimension, typename direction<Vector>::value_type> operator+(
		const direction<Vector>& u, const direction<Vector>& v) {
	typedef vector<vector_traits<Vector>::Dimension,
			typename direction<Vector>::value_type> result_type;
	result_type r;
	std::transform(begin(u), end(u), begin(v), begin(r),
			std::plus<typename result_type::value_type>());
	return r;
}

namespace inner {

template<typename Vector>
bool is_equal_impl(const direction<Vector>& u, const direction<Vector>& v,
		dimension<GK::GK_2D>) {
	return (u[GK::X] == v[GK::X] && u[GK::Y] == v[GK::Y]);
}

template<typename Vector>
bool is_equal_impl(const direction<Vector>& u, const direction<Vector>& v,
		dimension<GK::GK_3D>) {
	return (u[GK::X] == v[GK::X] && u[GK::Y] == v[GK::Y] && u[GK::Z] == v[GK::Z]);
}

}  // namespace inner

template<typename Vector>
bool operator==(const direction<Vector>& u, const direction<Vector>& v) {
	return inner::is_equal_impl(u, v,
			dimension<vector_traits<Vector>::Dimension>());
}

template<typename Vector>
bool operator!=(const direction<Vector>& u, const direction<Vector>& v) {
	return !(u == v);
}

template<typename T>
direction<geometry_traits<T>::Dimension> direction_of(const T& a);

template<typename Vector>
direction<vector_traits<Vector>::Dimension> normalize(const Vector& v) {
	return direction<vector_traits<Vector>::Dimension>(v);
}

template<typename Vector>
typename direction<Vector>::value_type norm(const direction<Vector>&) {
	return direction<Vector>::value_type(GK_FLOAT_ONE);
}

template<typename Vector>
Vector operator*(const typename vector_traits<Vector>::value_type& alpha,
		const direction<vector_traits<Vector>::Dimension>& u) {
//	return inner::multiply_scalar_direction(alpha, u,
//			dimension<direction<Vector>::Dimension>());
	Vector r;
	std::transform(begin(u), end(u), begin(r),
			std::bind2nd(
					multiplies<typename direction<Vector>::value_type,
							typename vector_traits<Vector>::value_type>(),
					alpha));
}

template<typename Vector>
Vector operator*(const direction<vector_traits<Vector>::Dimension>& u,
		const typename vector_traits<Vector>::value_type& alpha) {
	return operator*(alpha, u);
}

/**
 * @brief Basis in a vector space.
 * @author Takuya Makimoto
 * @date 2015/12/09
 */
template<size_t DimensionSize>
class basis {
public:
	static const size_t Dimension = DimensionSize;

public:
	basis() :
			X_() {
		this->initialize_();
	}

	basis(const basis& other) :
			X_() {

	}

	~basis() {

	}

	direction<DimensionSize> operator[](size_t n) const {
		return this->X_[n];
	}

private:
	direction<DimensionSize> X_[Dimension];

private:
	void initialize_() {
		for (size_t i = 0; i < Dimension; ++i) {
			this->X_[i][i] = GK_FLOAT_ONE;
		}
	}
};

}  // namespace gk

#endif /* GKDIRECTION_H_ */
