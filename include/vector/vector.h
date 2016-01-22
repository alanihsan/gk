/*
 * vector.h
 *
 *  Created on: 2016/01/22
 *      Author: tmakimoto
 */

#ifndef INCLUDE_VECTOR_VECTOR_H_
#define INCLUDE_VECTOR_VECTOR_H_

#include "../gkdef.h"

namespace gk {

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

template<size_t Dimension, typename T, typename X>
vector<Dimension, typename divides_result<T, X>::value_type> operator/(
		const vector<Dimension, T>& v, const X& alpha) {
	const T inv_alpha = alpha;
	return v * inv_alpha;
}

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

template<size_t Dimension, typename T>
vector<Dimension, typename direction<Dimension>::value_type> operator+(
		const direction<Dimension>& u, const direction<Dimension>& v) {
	typedef vector<Dimension, typename direction<Dimension>::value_type> result_type;
	result_type r;
	std::transform(begin(u), end(u), begin(v), begin(r),
			std::plus<typename result_type::value_type>());
	return r;
}

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

template<size_t DimensionSize, typename T>
struct vector_traits<vector<DimensionSize, T> > {
	typedef T value_type; ///< Type of elements in a vector.

	typedef value_type* iterator;
	typedef const value_type* const_iterator;
	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	static const size_t Dimension = DimensionSize; ///< A dimension size of a vector space.
	static const bool IsHomogeneous = false;

	typedef direction<Dimension> direction;
};

template<size_t DimensionSize>
struct vector_traits<direction<DimensionSize> > {
	typedef direction<DimensionSize> value_type; ///< Type of elements in a vector.

	typedef value_type* iterator;
	typedef const value_type* const_iterator;
	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	static const size_t Dimension = DimensionSize; ///< A dimension size of a vector space.
	static const bool IsHomogeneous = false;

	typedef direction<Dimension> direction;
};

} // namespace gk

#endif /* INCLUDE_VECTOR_VECTOR_H_ */
