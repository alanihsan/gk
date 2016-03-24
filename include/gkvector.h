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
#include "gkquaternion.h"

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
 * @tparam Vector A type of a vector in a vector space.
 *
 * @date 2016/02/24
 * @author Takuya Makimoto
 */
template<typename Vector>
struct vector_traits {
	typedef typename Vector::value_type value_type; ///< Type of elements in a vector.

	static const size_t Dimension = Vector::Dimension; ///< A dimension size of a vector space.
	static const bool IsHomogeneous = false;
};

template<typename Vector1, typename Vector2>
typename multiplies_result<typename vector_traits<Vector1>::value_type,
		typename vector_traits<Vector2>::value_type>::value_type dot(
		const Vector1& u, const Vector2& v) {

	typedef typename multiplies_result<
			typename vector_traits<Vector1>::value_type,
			typename vector_traits<Vector2>::value_type>::value_type result_type;

	const size_t Dimension = vector_traits<Vector1>::Dimension;

	result_type r = result_type(GK_FLOAT_ZERO);
	for (size_t i = 0; i < Dimension; ++i) {
		r += u[i] * v[i];
	}

	return r;
}

struct dot_function {
	template<typename Vector1, typename Vector2>
	typename multiplies_result<typename vector_traits<Vector1>::value_type,
			typename vector_traits<Vector2>::value_type>::value_type operator()(
			const Vector1& u, const Vector2& v) {
		typedef typename multiplies_result<
				typename vector_traits<Vector1>::value_type,
				typename vector_traits<Vector2>::value_type>::value_type result_type;

		const size_t Dimension = vector_traits<Vector1>::Dimension;

		result_type r = result_type(GK_FLOAT_ZERO);
		for (size_t i = 0; i < Dimension; ++i) {
			r += u[i] * v[i];
		}

		return r;
	}
};

template<typename Vector>
typename vector_traits<Vector>::value_type norm(const Vector& v) {
	return std::sqrt(dot(v, v));
}

namespace impl {

template<typename Vector1, typename Vector2, typename Result>
void cross_kernel(const Vector1&, const Vector2&, Result& r,
		dimension_tag<GK::GK_2D>) {
	typedef typename vector_traits<Result>::value_type value_type;
	r[GK::X] = value_type(GK_FLOAT_ZERO);
	r[GK::Y] = value_type(GK_FLOAT_ZERO);
}

template<typename Vector1, typename Vector2, typename Result>
void cross_kernel(const Vector1& u, const Vector2& v, Result& r,
		dimension_tag<GK::GK_3D>) {
	r[GK::X] = u[GK::Y] * v[GK::Z] - u[GK::Z] * v[GK::Y];
	r[GK::Y] = u[GK::Z] * v[GK::X] - u[GK::X] * v[GK::Z];
	r[GK::Z] = u[GK::X] * v[GK::Y] - u[GK::Y] * v[GK::X];
}

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
 * @brief Direction.
 *
 * @author Takuya Makimoto
 * @date 2016/01/25
 */
template<std::size_t DimensionSize>
class direction {
public:
	typedef gkfloat value_type;
	typedef const value_type* const_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	static const std::size_t Dimension = DimensionSize;
	static const std::size_t ElementSize = DimensionSize;

private:

//	/**
//	 * @brief Normalize a vector.
//	 * @param v
//	 * @param d
//	 */
//	template<typename Vector, size_t Dimension>
//	void Normalize_(const Vector& v, direction& d, dimension_tag<Dimension>) {
//		const typename divides_result<gkfloat,
//				typename vector_traits<Vector>::value_type>::value_type F =
//				gkfloat(
//				GK_FLOAT_ONE) / norm(v);
//		for (size_t i = 0; i < DimensionSize; ++i) {
//			d.x_[i] = F * v[i];
//		}
//	}
//
//	template<typename Vector>
//	void Normalize_(const Vector& v, direction& d, dimension_tag<GK::GK_2D>) {
//		const typename divides_result<gkfloat,
//				typename vector_traits<Vector>::value_type>::value_type F =
//				gkfloat(
//				GK_FLOAT_ONE) / norm(v);
//		d.x_[GK::X] = F * v[GK::X];
//		d.x_[GK::Y] = F * v[GK::Y];
//	}
//
//	template<typename Vector>
//	void Normalize_(const Vector& v, direction& d, dimension_tag<GK::GK_3D>) {
//		const typename divides_result<gkfloat,
//				typename vector_traits<Vector>::value_type>::value_type F =
//				gkfloat(
//				GK_FLOAT_ONE) / norm(v);
//		d.x_[GK::X] = F * v[GK::X];
//		d.x_[GK::Y] = F * v[GK::Y];
//		d.x_[GK::Z] = F * v[GK::Z];
//	}

	template<typename InputIterator, typename OutputIterator>
	void copy_(InputIterator first, OutputIterator result,
			std::input_iterator_tag) {
		InputIterator last = first;
		std::advance(last, DimensionSize);
		std::copy(first, last, result);
	}

	template<typename InputIterator, typename OutputIterator>
	void copy_(InputIterator first, OutputIterator result,
			std::random_access_iterator_tag) {
		std::copy(first, first + DimensionSize, result);
	}

	template<typename InputIterator, typename OutputIterator>
	void normalize_(InputIterator first, OutputIterator result) {
		typedef typename multiplies_result<
				typename std::iterator_traits<InputIterator>::value_type,
				typename std::iterator_traits<InputIterator>::value_type>::value_type T;

		InputIterator last = first;
		std::advance(last, DimensionSize);
		const typename divides_result<gkfloat,
				typename std::iterator_traits<InputIterator>::value_type>::value_type F =
				gkfloat(GK_FLOAT_ONE)
						/ std::sqrt(
								std::inner_product(first, last, first,
										T(GK_FLOAT_ZERO)));

//		std::transform(first,last,result,std::bind2nd(std::multiplies<gkfloat>))
	}

public:
	/**
	 * @brief Default constructor.
	 *
	 * An instance made by this constructor is implemented a zero vector.
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

	template<typename InputIterator>
	explicit direction(InputIterator x) :
			x_() {

	}

//	template<typename Vector>
//	explicit direction(const Vector& v) :
//			x_() {
//		Normalize_(v, *this, dimension_tag<DimensionSize>());
//	}

	~direction() {
	}

//	const_iterator begin() const {
//		return this->x_;
//	}
//
//	const_iterator end() const {
//		return this->x_ + Dimension;
//	}
//
//	const_reverse_iterator rbegin() const {
//		return std::reverse_iterator<const_iterator>(this->end());
//	}
//
//	const_reverse_iterator rend() const {
//		return std::reverse_iterator<const_iterator>(this->begin());
//	}

	const value_type* data() const {
		return this->x_;
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

private:
	value_type x_[ElementSize];
};

/**
 * @brief
 *
 * @author Takuya Makimoto
 * @date 2016/03/07
 */
template<size_t DimensionSize>
struct vector_traits<direction<DimensionSize> > {
	typedef typename direction<DimensionSize>::value_type value_type; ///< Type of elements in a vector.

	static const size_t Dimension = DimensionSize; ///< A dimension size of a vector space.
	static const bool IsHomogeneous = false;
};

template<size_t Dimension>
bool operator==(const direction<Dimension>& u, const direction<Dimension>& v) {
	return std::equal(u.begin(), u.end(), v.begin());
}

template<size_t Dimension>
bool operator!=(const direction<Dimension>& u, const direction<Dimension>& v) {
	return !(u == v);
}

template<typename CharT, typename Traits, size_t Dimension>
std::basic_ostream<CharT, Traits>& operator<<(
		std::basic_ostream<CharT, Traits>& os, const direction<Dimension>& v) {
	const std::streamsize n = os.width();
	for (typename direction<Dimension>::const_iterator p = v.begin();
			p != v.end(); ++p) {
		os.width(n);
		os << *p;
	}

	return os;
}

/**
 * @brief Outputs a direction of a geometry, @a a.
 *
 * @tparam T The type of the geometry.
 * @param a
 * @return
 */
template<typename T>
direction<geometry_traits<T>::Dimension> direction_of(const T& a);

/**
 * @brief Computes a direction of a vector.
 * @param v The vector to be computed.
 * @return
 */
template<typename Vector>
direction<vector_traits<Vector>::Dimension> normalize(const Vector& v) {
	return direction<vector_traits<Vector>::Dimension>(v);
}

/**
 * @brief Computes a norm of a direction. This function always returns @b1
 * because a direction is a unit vector.
 * @tparam Dimension Dimension of a vector space.
 * @param
 * @return Returns the magnitude of the direction, 1.
 */
template<size_t Dimension>
typename direction<Dimension>::value_type norm(const direction<Dimension>&) {
	return typename direction<Dimension>::value_type(GK_FLOAT_ONE);
}

namespace impl {

template<size_t Dimension, typename Vector>
direction<Dimension> gk_normal_direction(const Vector&, const Vector&,
		dimension_tag<Dimension>) {
	return direction<Dimension>();
}

template<typename Vector>
direction<GK::GK_3D> gk_normal_direction(const Vector& u, const Vector& v,
		dimension_tag<GK::GK_3D>) {
	const typename vector_traits<Vector>::value_type Unit(GK_FLOAT_ONE);

	Vector r;
	r[GK::X] = (u[GK::Y] * v[GK::Z] - u[GK::Z] * v[GK::Y]) / Unit;
	r[GK::Y] = (u[GK::Z] * v[GK::X] - u[GK::X] * v[GK::Z]) / Unit;
	r[GK::Z] = (u[GK::X] * v[GK::Y] - u[GK::Y] * v[GK::X]) / Unit;

	return direction<GK::GK_3D>(r);
}

} // namespace inner

/**
 * @brief Computes a normal vector made by 2 vectors @a u and @a v
 * in a 3D space.
 * @param u A vector for computing the normal vector.
 * @param v The other vector.
 * @return The normal vector to have been computed.
 */
template<typename Vector>
direction<vector_traits<Vector>::Dimension> normal_direction(const Vector& u,
		const Vector& v) {
	return impl::gk_normal_direction(u, v,
			dimension_tag<vector_traits<Vector>::Dimension>());
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
	typedef direction<DimensionSize> direction_type;

private:
	static direction_type Value_(size_t n) {
		gkfloat x[DimensionSize] = { gkfloat(GK_FLOAT_ZERO) };
		x[n] = GK_FLOAT_ONE;
		return direction_type(x);
	}

public:
	basis() {
	}

	~basis() {
	}

	direction_type operator[](size_t n) const {
		return Value_(n);
	}

private:
	basis(const basis&);
	basis& operator=(const basis&);
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
		return rotate_(v, dimension_tag<vector_traits<Vector>::Dimension>());
	}

private:
	Vector rotate_(const Vector& v, dimension_tag<GK::GK_2D>) const {
		const gkfloat sin = std::sin(this->angle);
		const gkfloat cos = std::cos(this->angle);

		Vector r;
		r[GK::X] = v[GK::X] * cos - v[GK::Y] * sin;
		r[GK::Y] = v[GK::X] * sin + v[GK::Y] * cos;
		return r;
	}

	Vector rotate_(const Vector& v, dimension_tag<GK::GK_3D>) const {
		const gkfloat theta = norm(this->angle);
		const direction<GK::GK_3D> axis(this->angle);

		const gkfloat sin = std::sin(0.5 * theta);
		const gkfloat cos = std::cos(0.5 * theta);

		const quaternion Q(sin * axis[quaternion::X], sin * axis[quaternion::Y],
				sin * axis[quaternion::Z], cos);

		const typename vector_traits<Vector>::value_type unit =
				typename vector_traits<Vector>::value_type(GK_FLOAT_ONE);

		const quaternion P(v / unit, quaternion::value_type(GK_FLOAT_ZERO));

//		const quaternion R = conj(Q) * P * Q;
		const quaternion R = Q * P * conj(Q);

		Vector r;
		r[GK::X] = R[quaternion::X] * unit;
		r[GK::Y] = R[quaternion::Y] * unit;
		r[GK::Z] = R[quaternion::Z] * unit;

		return r;
	}
};

}  // namespace gk

#ifdef GK_USING_EIGEN
#include "eigen/eigen_vector.h"
#endif

#endif /* GKVECTOR_H_ */
