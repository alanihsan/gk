/*
 * point.h
 *
 *  Created on: 2015/10/30
 *      Author: makitaku
 */

#ifndef PRIMITIVE_POINT_H_
#define PRIMITIVE_POINT_H_

namespace gk {

/**
 * @brief
 *
 * @author Takuya Makiomto
 * @date 2016/10/20
 */
template<typename T, std::size_t DimensionSize>
class point {
public:
	typedef T value_type;
	static const std::size_t Dimension = DimensionSize;

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

	const T& operator[](std::size_t n) const {
		return this->x_[n];
	}

	T& operator[](std::size_t n) {
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
	vector_type<T, DimensionSize>::type x_;

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
