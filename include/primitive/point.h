/*
 * point.h
 *
 *  Created on: 2015/10/30
 *      Author: makitaku
 */

#ifndef PRIMITIVE_POINT_H_
#define PRIMITIVE_POINT_H_

namespace gk {

/*template<typename Vector>*/
template<std::size_t DimensionSize, typename T>
class point {
public:
//	typedef Vector vector_type;
//	typedef typename vector_traits<Vector>::value_type value_type;
	typedef T value_type;
	static const size_t Dimension = DimensionSize; //vector_traits<Vector>::Dimension;

public:
	point() :
			x_() {
	}

	point(const point& other) :
			x_(other.x_) {
	}

	template<typename Vector>
	point(const Vector& v) :
			x_() {
		this->assign_(v);
	}

	~point() {
	}

//	vector_type position() const {
//		return this->x_;
//	}

	const T& operator[](std::size_t n) const {
		return this->x_[n];
	}

	T& operator[](size_t n) {
		return this->x_[n];
	}

	point& operator=(const point& q) {
		if (*q == this) {
			return *this;
		}

		this->x_ = q.x_;
		return *this;
	}

	template<typename Vector>
	point& operator=(const Vector& v) {
		this->assign_(v);
	}

private:
	T x_[DimensionSize];

private:
	template<typename Vector>
	void assign_(const Vector& v) {
		for (std::size_t i = 0; i < DimensionSize; ++i) {
			this->x_[i] = v[i];
		}
	}
};

}  // namespace gk

#endif /* PRIMITIVE_POINT_H_ */
