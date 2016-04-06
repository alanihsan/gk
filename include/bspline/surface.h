/*
 * surface.h
 *
 *  Created on: 2016/04/06
 *      Author: tmakimoto
 */

#ifndef INCLUDE_BSPLINE_SURFACE_H_
#define INCLUDE_BSPLINE_SURFACE_H_

#include "basis.h"

namespace gk {

template<typename Vector>
class network {
public:

	network() :
			column_size_(), Q_() {
	}

	network(const network& other) :
			column_size_(other.column_size_), Q_(other.Q_) {
	}

	network(std::size_t row_size, std::size_t column_size) :
			column_size_(column_size), Q_(row_size * column_size) {
	}

	template<typename InputIterator>
	network(InputIterator first, InputIterator last, std::size_t column_size) :
			column_size_(column_size), Q_(first, last) {
	}

	~network() {
	}

	std::size_t row_size() const {
		return this->row_size_();
	}

	std::size_t column_size() const {
		return this->column_size_;
	}

	const Vector& operator()(std::size_t row, std::size_t column) const {
		return this->Q_[this->element_index_(row, column)];
	}

	Vector& operator()(std::size_t row, std::size_t column) {
		return this->Q_[this->element_index_(row, column)];
	}

private:
	std::size_t column_size_;
	std::vector<Vector> Q_;

private:
	std::size_t row_size_() const {
		return this->Q_.size() - this->column_size_;
	}

	std::size_t element_index_(std::size_t row, std::size_t column) const {
		return column + row * this->column_size_;
	}
};

/**
 * @brief B-spline surface.
 *
 * @date 2016/04/06
 */
template<typename Vector, typename Parameter>
class bsurface {
public:
	typedef Vector vector_type;

private:
	knotvector<Parameter> S_;
	knotvector<Parameter> T_;
	network<Vector> Q_;
};

}  // namespace gk

#endif /* INCLUDE_BSPLINE_SURFACE_H_ */
