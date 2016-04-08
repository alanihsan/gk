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
			major_size_(), Q_() {
	}

	network(const network& other) :
			major_size_(other.major_size_), Q_(other.Q_) {
	}

	network(std::size_t major_size, std::size_t minor_size) :
			major_size_(major_size), Q_(major_size * minor_size) {
	}

	template<typename InputIterator>
	network(InputIterator first, InputIterator last, std::size_t column_size) :
			major_size_(column_size), Q_(first, last) {
	}

	~network() {
	}

	std::size_t major_size() const {
		return this->major_size_;
	}

	std::size_t minor_size() const {
		return this->minor_size_();
	}

	const Vector& operator()(std::size_t major, std::size_t minor) const {
		return this->Q_[this->element_index_(major, minor)];
	}

	Vector& operator()(std::size_t major, std::size_t minor) {
		return this->Q_[this->element_index_(major, minor)];
	}

private:
	std::size_t major_size_;
	std::vector<Vector> Q_;

private:
	std::size_t minor_size_() const {
		return this->Q_.size() - this->major_size_;
	}

	std::size_t element_index_(std::size_t major, std::size_t minor) const {
		return major + minor * this->major_size_;
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

public:
	bsurface() :
			S_(), T_(), Q_() {
	}

	bsurface(const bsurface& other) :
			S_(other.S_), T_(other.T_), Q_(other.Q_) {
	}

	~bsurface() {
	}

	std::size_t major_degree() const {
		return this->major_degree_();
	}

	std::size_t minor_degree() const {
		return this->minor_degree_();
	}

	Vector operator()(const Parameter& s, const Parameter& t) const {
		std::vector<gkfloat> M;
		basis_function(this->major_degree_(), this->S_.begin(), this->S_.end(),
				s, std::inserter(M, M.begin()));

		std::vector<gkfloat> N;
		basis_function(this->minor_degree_(), this->T_.begin(), this->T_.end(),
				t, std::inserter(N, N.begin()));

	}

	bsurface& operator=(const bsurface& rhs) {
		if (&rhs == this) {
			return *this;
		}

		this->S_ = rhs.S_;
		this->T_ = rhs.T_;
		this->Q_ = rhs.Q_;

		return *this;
	}

private:
	knotvector<Parameter> S_; ///< The knot vector in major order.
	knotvector<Parameter> T_; ///< The knot vector in minor order.
	network<Vector> Q_; ///< The control points.

private:
	std::size_t major_degree_() const {
		return bspline_degree(this->Q_.major_size(), this->S_.size());
	}

	std::size_t minor_degree_() const {
		return bspline_degree(this->Q_.minor_size(), this->T_.size());
	}
};

}  // namespace gk

#endif /* INCLUDE_BSPLINE_SURFACE_H_ */
