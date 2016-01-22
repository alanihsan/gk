/*
 * gkvector.h
 *
 *  Created on: 2015/07/08
 *      Author: makitaku
 */

#ifndef GKVECTOR_H_
#define GKVECTOR_H_

#include "gkdef.h"
#include "gkfunctional.h"
#include "gkgeometry.h"

#include <numeric>
#include <cmath>

#include <utility>
#include <ostream>
#include <iterator>
#include <algorithm>
#include <functional>

namespace gk {

/**
 * @brief Traits of a vector.
 */
template<typename Vector>
struct vector_traits {
	typedef typename Vector::value_type value_type; ///< Type of elements in a vector.

	typedef value_type* iterator;
	typedef const value_type* const_iterator;
	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	static const size_t Dimension = Vector::Dimension; ///< A dimension size of a vector space.
	static const bool IsHomogeneous = false;

	typedef Vector::direction_type direction_type;
};

template<typename Vector>
typename vector_traits<Vector>::const_iterator begin(const Vector& v) {
	return v.data();
}

template<typename Vector>
typename vector_traits<Vector>::iterator begin(Vector& v) {
	return v.data();
}

template<typename Vector>
typename vector_traits<Vector>::const_iterator end(const Vector& v) {
	return v.data() + vector_traits<Vector>::Dimension;
}

template<typename Vector>
typename vector_traits<Vector>::iterator end(Vector& v) {
	return v.data() + vector_traits<Vector>::Dimension;
}

template<typename Vector>
typename vector_traits<Vector>::const_reverse_iterator rbegin(const Vector& v) {
	return std::reverse_iterator<typename vector_traits<Vector>::const_iterator>(
			end(v));
}

template<typename Vector>
typename vector_traits<Vector>::reverse_iterator rbegin(Vector& v) {
	return std::reverse_iterator<typename vector_traits<Vector>::iterator>(
			end(v));
}

template<typename Vector>
typename vector_traits<Vector>::const_reverse_iterator rend(const Vector& v) {
	return std::reverse_iterator<typename vector_traits<Vector>::const_iterator>(
			begin(v));
}

template<typename Vector>
typename vector_traits<Vector>::reverse_iterator rend(Vector& v) {
	return std::reverse_iterator<typename vector_traits<Vector>::iterator>(
			begin(v));
}

template<typename Vector1, typename Vector2,
		typename Result = typename multiplies_result<
				typename vector_traits<Vector1>::value_type,
				typename vector_traits<Vector2>::value_type>::value_type>
struct dot {
	typedef Vector1 first_argument_type;
	typedef Vector2 second_argument_type;
	typedef Result result_type;

	Result operator()(const Vector1& u, const Vector2& v) const {
		return std::inner_product(begin(u), end(u), begin(v), Result());
	}
};

template<typename Vector>
typename vector_traits<Vector>::value_type norm(const Vector& v) {
	typedef typename vector_traits<Vector>::value_type value_type;
	const dot<Vector, Vector,
			typename multiplies_result<value_type, value_type>::value_type> dot;
	return std::sqrt(dot(v, v));
}

template<typename Vector1, typename Vector2 = Vector1, typename Result = Vector1>
struct cross {
	typedef Vector1 first_argument_type;
	typedef Vector2 second_argument_type;
	typedef Result result_type;

	Result operator()(const Vector1& u, const Vector2& v) const {
		Result r;
		r[GK::X] = u[GK::Y] * v[GK::Z] - u[GK::Z] * v[GK::Y];
		r[GK::Y] = u[GK::Z] * v[GK::X] - u[GK::X] * v[GK::Z];
		r[GK::Z] = u[GK::X] * v[GK::Y] - u[GK::Y] * v[GK::X];
		return r;
	}
};

/**
 * @brief
 * @tparam Vector
 * @tparam Angle
 *
 * @author Takuya Makimoto
 * @date 2016/01/20
 */
template<typename Vector, typename Angle>
struct rotate {
	typedef Vector argument_type;
	typedef Vector result_type;

	const Angle angle;

	rotate() :
			angle() {
	}

	rotate(const Angle& theta) :
			angle(theta) {
	}

	Vector operator()(const Vector& v) const {
		return rotate_(v, dimension<vector_traits<Vector>::Dimension>());
	}

private:
	Vector rotate_(const Vector& v, dimension<GK::GK_2D>) const {
		const gkfloat sin = std::sin(this->angle);
		const gkfloat cos = std::cos(this->angle);

		Vector r;
		r[GK::X] = v[GK::X] * cos - v[GK::Y] * sin;
		r[GK::Y] = v[GK::X] * sin + v[GK::Y] * cos;
		return r;
	}

	Vector rotate_(const Vector& v, dimension<GK::GK_3D>) const {

	}
};

/**
 * @brief Basis in a vector space.
 * @author Takuya Makimoto
 * @date 2015/12/09
 */
template<typename Vector>
class basis {
public:
	static const size_t Dimension = vector_traits<Vector>::Dimension;
	typedef typename vector_traits<Vector>::direction direction;

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

	direction operator[](size_t n) const {
		return this->X_[n];
	}

private:
	direction X_[Dimension];

private:
	void initialize_() {
		for (size_t i = 0; i < Dimension; ++i) {
			this->X_[i][i] = GK_FLOAT_ONE;
		}
	}
};

}  // namespace gk

#endif /* GKVECTOR_H_ */
