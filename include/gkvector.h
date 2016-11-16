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

template<typename T, std::size_t Dimension>
struct vector_dependence {
	typedef Eigen::Matrix<T, 1, Dimension, Eigen::RowMajor> vector_type;
};

/**
 * @brief Traits of a vector.
 * @tparam Vector The type of the vector in a vector space.
 *
 * @date 2016/02/24
 * @author Takuya Makimoto
 */
template<typename Vector>
struct vector_traits {
	typedef typename Vector::value_type value_type; ///< Type of elements in a vector.

	static const std::size_t Dimension = Vector::Dimension; ///< The dimension size of the vector space.
	static const bool IsHomogeneous = false; ///<

	typedef value_type* iterator;
	typedef const value_type* const_iterator;
	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	static const_iterator begin(const Vector& v) {
		return v.begin();
	}

	static iterator begin(Vector& v) {
		return v.begin();
	}

	static const_iterator end(const Vector& v) {
		return v.end();
	}

	static iterator end(Vector& v) {
		return v.end();
	}

	static const_reverse_iterator rbegin(const Vector& v) {
		return std::reverse_iterator<const_iterator>(v.end());
	}

	static reverse_iterator rbegin(Vector& v) {
		return std::reverse_iterator<iterator>(v.end());
	}

	static const_reverse_iterator rend(const Vector& v) {
		return std::reverse_iterator<const_iterator>(v.begin());
	}

	static reverse_iterator rend(Vector& v) {
		return std::reverse_iterator<iterator>(v.begin());
	}
};

template<std::size_t DimensionSize, typename T>
struct vector_traits<T[DimensionSize]> {
	typedef T Vector[DimensionSize];
	typedef T value_type;

	static const std::size_t Dimension = DimensionSize;
	static const bool IsHomogeneous = false;

	typedef T* iterator;
	typedef const T* const_iterator;
	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	static const_iterator begin(const Vector& v) {
		return v;
	}

	static iterator begin(Vector& v) {
		return v;
	}

	static const_iterator end(const Vector& v) {
		return v + DimensionSize;
	}

	static iterator end(Vector& v) {
		return v + DimensionSize;
	}

	static const_reverse_iterator rbegin(const Vector& v) {
		return std::reverse_iterator<const_iterator>(end(v));
	}

	static reverse_iterator rbegin(Vector& v) {
		return std::reverse_iterator<iterator>(end(v));
	}

	static const_reverse_iterator rend(const Vector& v) {
		return std::reverse_iterator<const_iterator>(begin(v));
	}

	static reverse_iterator rend(Vector& v) {
		return std::reverse_iterator<iterator>(begin(v));
	}
};

template<typename Vector1, typename Vector2>
typename multiplies_result<typename vector_traits<Vector1>::value_type,
		typename vector_traits<Vector2>::value_type>::value_type dot(
		const Vector1& u, const Vector2& v) {

#ifdef GK_DEBUG
#endif

	typedef typename multiplies_result<
			typename vector_traits<Vector1>::value_type,
			typename vector_traits<Vector2>::value_type>::value_type result_type;

	const std::size_t Dimension = vector_traits<Vector1>::Dimension;

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

//namespace impl {
//
//template<typename Vector1, typename Vector2, typename Result>
//void cross_kernel(const Vector1&, const Vector2&, Result& r,
//		dimension_tag<GK::GK_2D>) {
//	typedef typename vector_traits<Result>::value_type value_type;
//	r[GK::X] = value_type(GK_FLOAT_ZERO);
//	r[GK::Y] = value_type(GK_FLOAT_ZERO);
//}
//
//template<typename Vector1, typename Vector2, typename Result>
//void cross_kernel(const Vector1& u, const Vector2& v, Result& r,
//		dimension_tag<GK::GK_3D>) {
//	r[GK::X] = u[GK::Y] * v[GK::Z] - u[GK::Z] * v[GK::Y];
//	r[GK::Y] = u[GK::Z] * v[GK::X] - u[GK::X] * v[GK::Z];
//	r[GK::Z] = u[GK::X] * v[GK::Y] - u[GK::Y] * v[GK::X];
//}
//
//}

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
template<std::size_t Dimension>
class direction {
public:
	typedef float_type value_type;
	typedef const value_type* const_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

//	static const std::size_t Dimension = DimensionSize;
//	static const std::size_t ElementSize = DimensionSize;

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
		std::copy(other.x_, other.x_ + Dimension, this->x_);
	}

//	template<typename InputIterator>
//	explicit direction(InputIterator first) :
//			x_() {
//		typedef typename std::iterator_traits<InputIterator>::value_type L_t;
//		typedef typename multiplies_result<L_t, L_t>::value_type L2_t;
//		typedef typename divides_result<gkfloat, L_t>::value_type InvL_t;
//
//		InputIterator last = first;
//		std::advance(last, DimensionSize);
//		const L2_t L2 = std::inner_product(first, last, first,
//				L2_t(GK_FLOAT_ZERO));
//		const InvL_t F = gkfloat(GK_FLOAT_ONE) / std::sqrt(L2);
//
//		std::transform(first, last, this->x_,
//				std::bind2nd(multiplies<L_t, InvL_t>(), F));
//	}

	template<typename Vector>
	direction(const Vector& v) :
			x_() {

	}

	template<typename Vector>
	direction(const Vector& start, const Vector& end) :
			x_() {

		typedef vector_traits<Vector> vtraits;
		typedef typename vtraits::value_type L_t;
		typedef typename multiplies_result<L_t, L_t>::value_type L2_t;
		typedef typename divides_result<value_type, L_t>::value_type InvL_t;

		Vector v = end - start;
		typename vtraits::iterator first = vtraits::begin(v);
		typename vtraits::iterator last = vtraits::end(v);
		const L2_t L2 = std::inner_product(first, last, first,
				L2_t(GK_FLOAT_ZERO));
		const InvL_t F = value_type(GK_FLOAT_ONE) / std::sqrt(L2);

		std::transform(first, last, this->x_,
				std::bind2nd(multiplies<L_t, InvL_t>(), F));
	}

	~direction() {
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

	value_type operator[](const size_t n) const {
		return this->x_[n];
	}

	direction& operator=(const direction& u) {
		if (&u == this) {
			return *this;
		}

		std::copy(u.x_, u.x_ + Dimension, this->x_);
		return *this;
	}

private:
	value_type x_[Dimension];
};

/**
 * @brief
 *
 * @date 2016/03/07
 */
template<size_t DimensionSize>
struct vector_traits<direction<DimensionSize> > {
	typedef typename direction<DimensionSize>::value_type value_type; ///< Type of elements in a vector.

	static const size_t Dimension = DimensionSize; ///< A dimension size of a vector space.
	static const bool IsHomogeneous = false;

	typedef typename direction<DimensionSize>::const_iterator iterator;
	typedef typename direction<DimensionSize>::const_iterator const_iterator;
	typedef typename direction<DimensionSize>::const_reverse_iterator reverse_iterator;
	typedef typename direction<DimensionSize>::const_reverse_iterator const_reverse_iterator;

	static const_iterator begin(const direction<DimensionSize>& d) {
		return d.begin();
	}

	static iterator begin(direction<DimensionSize>& d) {
		return d.begin();
	}

	static const_iterator end(const direction<DimensionSize>& d) {
		return d.end();
	}

	static iterator end(direction<DimensionSize>& d) {
		return d.end();
	}

	static const_reverse_iterator rbegin(const direction<DimensionSize>& d) {
		return d.rbegin();
	}

	static reverse_iterator rbegin(direction<DimensionSize>& d) {
		return d.rbegin();
	}

	static const_reverse_iterator rend(const direction<DimensionSize>& d) {
		return d.rend();
	}

	static reverse_iterator rend(direction<DimensionSize>& d) {
		return d.rend();
	}
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
 *
 * @f[
 * \mathbf{d} = \frac{\mathbf{v}}{|\mathbf{v}|}
 * @f]
 *
 * @param v The vector to be computed.
 * @return
 */
template<typename Vector>
direction<vector_traits<Vector>::Dimension> normalize(const Vector& v) {
	return direction<vector_traits<Vector>::Dimension>(
			vector_traits<Vector>::begin(v));
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
	static direction_type Value_(std::size_t n) {
		float_type x[DimensionSize] = { float_type(GK_FLOAT_ZERO) };
		x[n] = GK_FLOAT_ONE;
		return direction_type(x);
	}

public:
	basis() {
	}

	~basis() {
	}

	direction_type operator[](std::size_t n) const {
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
		const float_type sin = std::sin(this->angle);
		const float_type cos = std::cos(this->angle);

		Vector r;
		r[GK::X] = v[GK::X] * cos - v[GK::Y] * sin;
		r[GK::Y] = v[GK::X] * sin + v[GK::Y] * cos;
		return r;
	}

	Vector rotate_(const Vector& v, dimension_tag<GK::GK_3D>) const {
		const float_type theta = norm(this->angle);
		const direction<GK::GK_3D> axis(this->angle);

		const float_type sin = std::sin(0.5 * theta);
		const float_type cos = std::cos(0.5 * theta);

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

#if defined(GK_EIGEN_ROWVECTOR) || defined(GK_EIGEN_COLUMNVECTOR)
#include "eigen/eigen_vector.h"
#endif

#endif /* GKVECTOR_H_ */
