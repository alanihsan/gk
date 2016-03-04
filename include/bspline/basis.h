/*
 * basis.h
 *
 *  Created on: 2015/06/12
 *      Author: makitaku
 */

#ifndef BSPLINE_BASIS_H_
#define BSPLINE_BASIS_H_

#include <gkdef.h>
#include <vector>
#include <algorithm>

namespace gk {

size_t bspline_degree(size_t control_points_size, size_t knot_vector_size) {
	return knot_vector_size - control_points_size - 1;
}

size_t bspline_control_points_size(size_t degree, size_t knot_vector_size) {
	return knot_vector_size - degree - 1;
}

size_t bspline_knot_vector_size(size_t degree, size_t control_points_size) {
	return control_points_size + degree + 1;
}

/**
 * @brief Knot vector.
 *
 * @author Takuya Makimoto
 * @date 2016/01/07
 */
template<typename T>
class knotvector {
public:
	typedef T value_type;
	typedef std::vector<T> container_type;

	typedef typename container_type::const_reference const_reference;
	typedef typename container_type::const_iterator const_iterator;
	typedef typename container_type::const_reverse_iterator const_reverse_iterator;

public:
	knotvector() :
			T_() {
	}

	knotvector(const knotvector& other) :
			T_(other.T_) {
	}

	knotvector(gksize size) :
			T_(size) {
	}

	template<typename InputIterator>
	knotvector(InputIterator first, InputIterator last) :
			T_(first, last) {
		std::stable_sort(this->T_.begin(), this->T_.end());
	}

	~knotvector() {
	}

	bool empty() const {
		return this->T_.empty();
	}

	gksize size() const {
		return this->T_.size();
	}

	const_reference front() const {
		return this->T_.front();
	}

	const_reference back() const {
		return this->T_.back();
	}

	const_iterator begin() const {
		return this->T_.begin();
	}

	const_iterator end() const {
		return this->T_.end();
	}

	const_reverse_iterator rbegin() const {
		return this->T_.rbegin();
	}

	const_reverse_iterator rend() const {
		return this->T_.rend();
	}

	const_iterator insert(const value_type& t) {
		typename container_type::iterator p = std::upper_bound(this->T_.begin(),
				this->T_.end(), t);
		return this->T_.insert(p, t);
	}

	void insert(gkfloat t, gksize multiplicity) {
		typename container_type::iterator p = std::upper_bound(this->T_.begin(),
				this->T_.end(), t);
		this->T_.insert(p, multiplicity, t);

	}

	template<class InputIterator>
	void insert(InputIterator first, InputIterator last) {
		this->T_.insert(this->T_.end(), first, last);
		std::stable_sort(this->T_.begin(), this->T_.end());
	}

	/**
	 * @brief Erases the element at @a position.
	 * @param position
	 * @return
	 */
	const_iterator erase(const_iterator position) {
		return this->T_.erase(this->T_.begin() + (position - this->T_.begin()));
	}

	/**
	 * @brief Erases the elements between @a first and the previous position of
	 * @a last, [first, last).
	 * @param first
	 * @param last
	 * @return
	 */
	const_iterator erase(const_iterator first, const_iterator last) {
		return this->T_.erase(this->T_.begin() + (first - this->T_.begin()),
				this->T_.begin() + (last - this->T_.begin()));
	}

	value_type operator[](gksize n) const {
		return this->T_[n];
	}

	knotvector& operator=(const knotvector& rhs) {
		if (&rhs == this) {
			return *this;
		}

		this->T_ = rhs.T_;
		return *this;
	}

private:
	container_type T_;
};

template<typename Parameter>
size_t knot_segment(const knotvector<Parameter>& T, std::size_t degree,
		const Parameter& t) {
	const std::size_t order = degree + 1;
	const std::size_t n = T.size() - order;

	return std::distance(T.begin(),
			std::upper_bound(T.begin() + order, T.begin() + n, t)) - 1;
}

/**
 * @brief Compute basis function values at a parameter @a t.
 * @param degree Degree of the basis function.
 * @param T Object of knot vector.
 * @param t Parameter.
 * @param N Beginning iterator of the basis function.
 * @param dorder Derivative order.
 *
 * @return The last iterator of the basis function.
 */
template<typename InputRandomAccessIterator, typename Parameter,
		typename OutputRandomAccessIterator>
OutputRandomAccessIterator basis_function(size_t degree,
		InputRandomAccessIterator first, InputRandomAccessIterator last,
		const Parameter& t, OutputRandomAccessIterator N, size_t dorder = 0) {

	const Parameter Zero = Parameter(GK_FLOAT_ZERO);
	const Parameter One = Parameter(GK_FLOAT_ONE);

	const size_t order = degree + 1;
	const size_t knot_size = std::distance(first, last);
	const size_t control_size = knot_size - order;

	const size_t segment = std::distance(first,
			std::upper_bound(first + order, first + control_size, t)) - 1;
	InputRandomAccessIterator T = first;

	/*
	 * Initialization.
	 */
	std::fill_n(N, control_size, Zero);

	N[segment] = One;

	if (degree == 0) {
		return ++N;
	}

	typedef std::vector<Parameter> coeff_t;
	for (size_t k = 2; k <= order - dorder; ++k) {

		coeff_t alpha(control_size);

		for (size_t j = segment - k + 2; j <= segment; j++) {
			const Parameter dt = T[j + k - 1] - T[j];
			alpha[j] = (dt == Zero) ? Zero : (t - T[j]) / dt;
		}

		N[segment - k + 1] = (One - alpha[segment - k + 2])
				* N[segment - k + 2];
		for (size_t j = segment - k + 2; j < segment; j++) {
			N[j] = alpha[j] * N[j] + (One - alpha[j + 1]) * N[j + 1];
		}
		N[segment] = alpha[segment] * N[segment];
	}

//	return N + control_size;

	if (dorder == 0) {
		return N + control_size;
	} else {

		for (size_t k = order - dorder + 1; k <= order; ++k) {
			const size_t p = k - 1; // degree

			coeff_t beta(control_size);

			for (size_t j = segment - k + 2; j <= segment; ++j) {
				const Parameter denominator = T[j + p] - T[j];
				beta[j] = (denominator == Zero) ? Zero : One / denominator;
			}

			N[segment - k + 1] = -p * beta[segment - k + 2]
					* N[segment - k + 2];
			for (gksize j = segment - k + 2; j < segment; ++j) {
				N[j] = p * (beta[j] * N[j] - beta[j + 1] * N[j + 1]);
			}
			N[segment] = p * beta[segment] * N[segment];
		}

		return N + control_size;
	}
}

//template<typename Vector>
//class bspline_basis {
//public:
//	typedef Vector vector_type;
//	typedef typename vector_traits<Vector>::value_type value_type;
//
//public:
//	bspline_basis() :
//			derivatives_() {
//	}
//
//	bspline_basis(const bspline_basis& other) :
//			derivatives_(other.derivatives_) {
//	}
//
//	~bspline_basis() {
//	}
//
//	vector_type operator()(const knotvector<value_type>& T, std::size_t degree,
//			const value_type& t) const {
//		const value_type zero(GK_FLOAT_ZERO);
//		const value_type one(GK_FLOAT_ZERO);
//
//		const std::size_t order = degree + 1;
//		const std::size_t n = T.size() - order; // number of control points.
//		vector_type N(n);
//
//		const std::size_t segment = knot_segment(T, degree, t);
//
//		N[segment] = one; // degree = 0.
//
//		if (degree == 0) {
//			return N;
//		}
//
//		const std::size_t D = degree - this->derivatives_;
//		for (std::size_t p = 1; p <= D; ++p) {
//			vector_type alpha(n);
//			for (std::size_t j = segment - p + 1; j <= segment; ++j) {
//				const value_type denominator = T[j + p] - T[j];
//				alpha[j] =
//						(denominator == zero) ? zero : (t - T[j]) / denominator;
//			}
//
//			N[segment - p] = (one - alpha[segment - p + 1])
//					* N[segment - p + 1];
//			for (std::size_t j = segment - p + 1; j < segment; ++j) {
//				N[j] = alpha[j] * N[j] + (one - alpha[j + 1]) * N[j + 1];
//			}
//			N[segment] = alpha[segment] * N[segment];
//		}
//
//		if (this->derivatives_ == 0) {
//			return N;
//		}
//
//		for (std::size_t p = D + 1; p <= degree; ++p) {
//			vector_type beta(n);
//
//			for (std::size_t j = segment - p + 1; j <= segment; ++j) {
//				const value_type denominator = T[j + p] - T[j];
//				beta[j] =
//						(denominator == zero) ?
//								zero :
//								static_cast<typename scalar_traits<value_type>::float_type>(p)
//										/ denominator;
//			}
//
//			N[segment - p] = -beta[segment - p + 1] * N[segment - p + 1];
//			for (std::size_t j = segment - p + 1; j < segment; ++j) {
//				N[j] = beta[j] * N[j] - beta[j + 1] * N[j + 1];
//			}
//			N[segment] = beta[segment] * N[segment];
//		}
//
//		return N;
//	}
//
//private:
//	const std::size_t derivatives_;
//};
//
//gksize search_segment(const gkknotvector& T, gksize order, gkfloat t) {
//	const gksize control_size = T.size() - order;
//
//	return std::distance(T.begin(),
//			std::upper_bound(T.begin() + order, T.begin() + control_size, t));
//}
//
//template<gksize DerivativeOrder> class basis;
//
//template<>
//class basis<0> {
//public:
//	typedef Eigen::Matrix<gkfloat, 1, Eigen::Dynamic> vector;
//public:
//	basis(const basis& other) :
//			T_(other.T_), order_(other.order_) {
//
//	}
//
//	explicit basis(const gkknotvector& T, gksize order) :
//			T_(T), order_(order) {
//
//	}
//
//	~basis() {
//
//	}
//
//	vector operator()(gkfloat t) const {
//		return operator()(this->order_, t);
//	}
//
//	vector operator()(gksize order, gkfloat t) const {
//		const gksize control_size = this->T_.size() - this->order_;
//		typedef vector dest_type;
//		dest_type N(control_size);
//
////		const gksize degree = this->order_ - 1;
//
////			gkknotvector::const_iterator p = basis::T_.segment(this->t_);
//		const typename gkknotvector::const_iterator first = this->T_.begin()
//				+ this->order_;
//		const typename gkknotvector::const_iterator end = this->T_.begin()
//				+ control_size;
//		const gksize segment = std::distance(this->T_.begin(),
//				std::upper_bound(first, end, t)) - 1;
//
//		N[segment] = GK::Unit;
//
//		for (gksize k = 2; k <= order; k++) {
//
//			typedef dest_type coeff_type;
//			coeff_type alpha(control_size);
//			//				gksize j = order - k;
//
//			for (gksize j = segment - k + 1; j <= segment; j++) {
//				const gkfloat denominator = this->T_[j + k - 1] - this->T_[j];
//				alpha[j] =
//						(denominator == GK::Zero) ?
//								GK::Zero : (t - this->T_[j]) / denominator;
//			}
//
//			N[segment - k + 1] = (GK::Unit - alpha[segment - k + 2])
//					* N[segment - k + 2];
//			for (gksize j = segment - k + 2; j < segment; j++) {
//				N[j] = alpha[j] * N[j] + (GK::Unit - alpha[j + 1]) * N[j + 1];
//			}
//			N[segment] = alpha[segment] * N[segment];
//		}
//
//		return N;
//	}
//
//private:
//	const gkknotvector& T_;
//	const gksize order_;
//
//private:
//	basis();
//	basis& operator=(const basis&);
//
////	void calcurate_coefficient_(gkfloat t, gksize order, gksize segment,
////			Eigen::Matrix<gkfloat, 1, Eigen::Dynamic>& N) {
////		const gksize degree = order - 1;
////
////		Eigen::Matrix<gkfloat, 1, Eigen::Dynamic> alpha(N.cols());
////
////		for (gksize j = segment - degree; j <= segment; j++) {
////			const gkfloat denominator = this->T_[j + degree] - this->T_[j];
////			alpha[j] =
////					(denominator == GK::Zero) ?
////							GK::Zero : (t - this->T_[j]) / denominator;
////		}
////
////		N[segment - degree] = (GK::Unit - alpha[segment - degree + 1])
////				* N[segment - degree + 1];
////		for (gksize j = segment - degree + 1; j < segment; j++) {
////			N[j] = alpha[j] * N[j] + (GK::Unit - alpha[j + 1] * N[j + 1]);
////		}
////		N[segment] *= alpha[segment];
////	}
//};
//
//template<gksize Derivative>
//class basis {
//public:
//	basis(const basis& other) :
//			T_(other.T_) {
//
//	}
//
//	basis(const gkknotvector& T, gksize order) :
//			T_(T), order_(order) {
//
//	}
//
//	~basis() {
//
//	}
//
//	Eigen::Matrix<gkfloat, 1, Eigen::Dynamic> operator()(gkfloat t) const {
//		return this->operator ()(this->order_, t);
//	}
//
//	Eigen::Matrix<gkfloat, 1, Eigen::Dynamic> operator()(gksize order,
//			gkfloat t) const {
//		typedef Eigen::Matrix<gkfloat, 1, Eigen::Dynamic> dest_type;
//
//		basis<0> N(this->T_, this->order_);
//
//		if (this->order_ - Derivative >= order) {
//			return N(order, t);
//		}
//
//		dest_type dN = N(this->order_ - Derivative, t);
//
//		const gksize control_size = dN.cols();
//
//		const gksize segment = std::distance(this->T_.begin(),
//				std::upper_bound(this->T_.begin() + this->order_,
//						this->T_.begin() + control_size, t)) - 1;
//
////		const gksize degree = this->order_ - 1;
//
//		typedef dest_type coeff_type;
//		coeff_type beta(control_size);
//
////		dest_type dN(control_size);
//
//		for (gksize k = this->order_ - Derivative + 1; k <= order; k++) {
//			const gksize p = k - 1; // degree
//
//			for (gksize j = segment - k + 2; j <= segment; j++) {
//				const gkfloat denominator = this->T_[j + p] - this->T_[j];
//				beta[j] =
//						(denominator == GK::Zero) ?
//								GK::Zero :
//								static_cast<gkfloat>(p) / denominator;
//			}
//
//			dN[segment - k + 1] = -beta[segment - k + 2] * dN[segment - k + 2];
//			for (gksize j = segment - k + 2; j < segment; j++) {
//				dN[j] = beta[j] * dN[j] - beta[j + 1] * dN[j + 1];
//			}
//			dN[segment] = beta[segment] * dN[segment];
//		}
//
//		return dN;
//	}
//
//private:
//	const gkknotvector& T_;
//	const gksize order_;
//
//private:
//	basis();
//	basis& operator=(const basis&);
//
//};

}// namespace gk

#endif /* BSPLINE_BASIS_H_ */
