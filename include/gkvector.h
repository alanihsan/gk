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
 * @brief Direction.
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
		std::transform(gk::begin(v), gk::end(v), d.x_,
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

//template<typename Vector>
//vector<vector_traits<Vector>::Dimension, typename direction<Vector>::value_type> operator+(
//		const direction<Vector>& u, const direction<Vector>& v) {
//	typedef vector<vector_traits<Vector>::Dimension,
//			typename direction<Vector>::value_type> result_type;
//	result_type r;
//	std::transform(begin(u), end(u), begin(v), begin(r),
//			std::plus<typename result_type::value_type>());
//	return r;
//}

template<size_t Dimension>
bool operator==(const direction<Dimension>& u, const direction<Dimension>& v) {
	return std::equal(u.begin(), u.end(), v.begin());
}

template<size_t Dimension>
bool operator!=(const direction<Dimension>& u, const direction<Dimension>& v) {
	return !(u == v);
}

template<typename T>
direction<geometry_traits<T>::Dimension> direction_of(const T& a);

template<typename Vector>
direction<vector_traits<Vector>::Dimension> normalize(const Vector& v) {
	return direction<vector_traits<Vector>::Dimension>(v);
}

template<size_t Dimension>
typename direction<Dimension>::value_type norm(const direction<Dimension>&) {
	return direction<Dimension>::value_type(GK_FLOAT_ONE);
}

//template<typename T>
//struct vector_from {
//	typedef void value_type;
//};

//template<typename X>
//typename vector_from<X>::value_type operator*(const X& alpha,
//		const direction<Dimension>& u) {
////	return inner::multiply_scalar_direction(alpha, u,
////			dimension<direction<Vector>::Dimension>());
//	typedef typename vector_from<X>::value_type vector;
//	vector r;
//	std::transform(begin(u), end(u), begin(r),
//			std::bind2nd(
//					multiplies<
//							typename direction<vector_traits<vector>::Dimension>::value_type,
//							typename vector_traits<vector>::value_type>(),
//					alpha));
//}
//
//template<typename Vector>
//Vector operator*(const direction<vector_traits<Vector>::Dimension>& u,
//		const typename vector_traits<Vector>::value_type& alpha) {
//	return operator*(alpha, u);
//}

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

/*
 * Declaration of vector_base.
 */
template<typename Derived>
class vector_base;

/**
 * @brief Vector. This is the class specified template class
 * vector for 3D.
 * @f[
 * \mathbf{x} =
 * \left[
 * \begin{array}{ccc}
 * x & y & z
 * \end{array}
 * \right]
 * @f]
 *
 * @author Takuya Makimoto
 * @date 2015
 */
template<size_t DimensionSize, typename T,
		template<size_t, typename > class Derived>
class vector_base<Derived<DimensionSize, T> > {
public:
	static const size_t Dimension = DimensionSize;

	typedef T value_type;
	typedef T* iterator;
	typedef const T* const_iterator;
	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

public:
	~vector_base() {

	}

	const_iterator begin() const {
		return this->x_;
	}

	iterator begin() {
		return this->x_;
	}

	const_iterator end() const {
		return this->x_ + Dimension;
	}

	iterator end() {
		return this->x_ + Dimension;
	}

	const_reverse_iterator rbegin() const {
		return const_reverse_iterator(this->x_ + Dimension);
	}

	reverse_iterator rbegin() {
		return reverse_iterator(this->x_ + Dimension);
	}

	const_reverse_iterator rend() const {
		return this->x_;
	}

	reverse_iterator rend() {
		return this->x_;
	}

	const T& operator[](size_t n) const {
		return this->x_[n];
	}

	T& operator[](size_t n) {
		return this->x_[n];
	}

	Derived<DimensionSize, T>& operator=(const vector_base& v) {
		if (&v == this) {
			return static_cast<Derived<DimensionSize, T>&>(*this);
		}

		std::copy(v.x_, v.x_ + Dimension, this->x_);
		return static_cast<Derived<DimensionSize, T>&>(*this);
	}

	Derived<DimensionSize, T>& operator+=(const vector_base& v) {
		std::transform(v.x_, v.x_ + DimensionSize, this->x_, this->x_,
				std::plus<T>());
		return static_cast<Derived<DimensionSize, T>&>(*this);
	}

	Derived<DimensionSize, T>& operator-=(const vector_base& v) {
		std::transform(v.x_, v.x_ + DimensionSize, this->x_, this->x_,
				std::minus<T>());
		return static_cast<Derived<DimensionSize, T>&>(*this);
	}

	Derived<DimensionSize, T>& operator*=(const T& alpha) {
		std::transform(this->x_, this->x_ + DimensionSize, this->x_,
				std::bind2nd(std::multiplies<T>(), alpha));
		return static_cast<Derived<DimensionSize, T>&>(*this);
	}

	Derived<DimensionSize, T>& operator/=(const T& alpha) {
		const T inv_alpha = T(GK_FLOAT_ONE) / alpha;
		return operator*=(inv_alpha);
	}

	const T* data() const {
		return this->x_;
	}

	T* data() {
		return this->x_;
	}

protected:
	vector_base() :
			x_() {

	}

	vector_base(const vector_base& other) :
			x_() {
		std::copy(other.x_, other.x_ + Dimension, this->x_);
	}

private:
	T x_[DimensionSize];
};

template<size_t Dimension, typename T>
class vector: public vector_base<vector<Dimension, T> > {
public:
	typedef vector_base<vector> base;

public:
	vector() :
			base() {

	}

	vector(const vector& other) :
			base(other) {

	}

	~vector() {

	}
};

template<typename T>
class vector<GK::GK_2D, T> : public vector_base<vector<GK::GK_2D, T> > {
public:
	typedef vector_base<vector> base;
public:
	vector() :
			base() {

	}

	vector(const vector& other) :
			base(other) {

	}

	vector(const T& x, const T& y) :
			base() {
		base::operator[](GK::X) = x;
		base::operator[](GK::Y) = y;
	}

	~vector() {

	}
};

template<typename T>
class vector<GK::GK_3D, T> : public vector_base<vector<GK::GK_3D, T> > {
public:
	typedef vector_base<vector> base;

public:
	vector() :
			base() {

	}

	vector(const vector& other) :
			base(other) {

	}

	vector(const T& x, const T& y, const T& z) :
			base() {
		base::operator[](GK::X) = x;
		base::operator[](GK::Y) = y;
		base::operator[](GK::Z) = z;
	}

	~vector() {

	}
};

template<typename CharT, typename Traits, size_t Dimension, typename T>
std::basic_ostream<CharT, Traits>& operator<<(
		std::basic_ostream<CharT, Traits>& os, const vector<Dimension, T>& v) {
	const std::streamsize n = os.width();
	for (typename vector<Dimension, T>::const_iterator p = v.begin();
			p != v.end(); ++p) {
		os.width(n);
		os << *p;
	}

	return os;
}

template<size_t Dimension, typename T>
vector<Dimension, T> operator+(const vector<Dimension, T>& u,
		const vector<Dimension, T>& v) {
	vector<Dimension, T> w = u;
	w += v;
	return w;
}

template<size_t Dimension, typename T>
vector<Dimension, T> operator-(const vector<Dimension, T>& u,
		const vector<Dimension, T>& v) {
	vector<Dimension, T> w = u;
	w -= v;
	return w;
}

template<size_t Dimension, typename T, typename X>
vector<Dimension, typename multiplies_result<T, X>::value_type> operator*(
		const vector<Dimension, T>& v, const X& alpha) {
	typedef typename multiplies_result<T, X>::value_type result_value;
	typedef vector<Dimension, typename multiplies_result<T, X>::value_type> result_vector;

	result_vector w;

	std::transform(begin(v), end(v), begin(w),
			std::bind2nd(multiplies<T, X, result_value>(), alpha));
	return w;
}

template<size_t Dimension, typename T, typename X>
vector<Dimension, typename multiplies_result<T, X>::value_type> operator*(
		const X& alpha, const vector<Dimension, T>& v) {
	return v * alpha;
}

template<size_t Dimension, typename T>
vector<Dimension, T> operator*(const T& alpha, const direction<Dimension>& u) {
	vector<Dimension, T> v;
	std::transform(u.begin(), u.end(), v.begin(),
			std::bind2nd(
					multiplies<T, typename direction<Dimension>::value_type>(),
					alpha));

	return v;
}

template<size_t Dimension, typename T>
vector<Dimension, T> operator*(const direction<Dimension>& u, const T& alpha) {
	return alpha * u;
}

template<size_t Dimension, typename T, typename X>
vector<Dimension, typename divides_result<T, X>::value_type> operator/(
		const vector<Dimension, T>& v, const X& alpha) {
	const T inv_alpha = alpha;
	return v * inv_alpha;
}

}  // namespace gk

#endif /* GKVECTOR_H_ */
