/*
 * point.h
 *
 *  Created on: 2015/10/30
 *      Author: makitaku
 */

#ifndef PRIMITIVE_POINT_H_
#define PRIMITIVE_POINT_H_

namespace gk {

template<typename Vector>
class point/*: geometry<point_tag, Vector,
 typename vector_traits<Vector>::value_type>*/{
public:
	typedef Vector vector_type;
	typedef typename vector_traits<Vector>::value_type value_type;

	static const size_t Dimension = vector_traits<Vector>::Dimension;

public:
	point() :
			x_() {

	}

	point(const point& other) :
			x_(other.x_) {

	}

	point(const vector_type& v) :
			x_(v) {

	}

	~point() {

	}

	vector_type position() const {
		return this->x_;
	}

	const value_type& operator[](size_t n) const {
		return this->x_[n];
	}

	value_type& operator[](size_t n) {
		return this->x_[n];
	}

	point& operator=(const point& q) {
		if (*q == this) {
			return *this;
		}

		this->x_ = q.x_;
		return *this;
	}

private:
	vector_type x_;
};

}  // namespace gk

#endif /* PRIMITIVE_POINT_H_ */
