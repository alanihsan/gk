/*
 * bspline.h
 *
 *  Created on: 2014/07/16
 *      Author: makitaku
 */

#ifndef BSPLINE_BSPLINE_H_
#define BSPLINE_BSPLINE_H_

#include <vector>
#include <utility>
#include <algorithm>
#include <numeric>

#include <gkvector.h>
#include <gkcurve.h>
#include <gkaabb.h>
#include "../primitive/line.h"
#include "basis.h"

namespace gk {

/**
 * @brief B-spline (Basis spline).
 *
 * @tparam Vector Type of a control point.
 *
 * @author Takuya Makimoto
 * @date 2015
 */
template<typename Vector, typename Parameter = gkfloat>
class bspline: public curve<free_curve_tag, Vector, Parameter> {
public:
	typedef curve<free_curve_tag, Vector,
			typename vector_traits<Vector>::value_type> base_type;
	typedef typename base_type::vector_type vector_type;
	typedef typename base_type::parameter parameter;
	typedef knotvector<parameter> knotvector_type;
	typedef std::vector<Vector> control_points;

private:
	static size_t compute_degree_(size_t knotvector_size,
			size_t controls_size) {
		return knotvector_size - controls_size - 1;
	}

public:

	bspline() :
			T_(), Q_() {
	}

	bspline(const bspline& other) :
			T_(other.T_), Q_(other.Q_) {

	}

	bspline(const knotvector_type& T, const control_points& Q) :
			T_(T), Q_(Q) {
	}

//	template<typename InputIterator>
//	bspline(const knotvector<parameter>& T, InputIterator first,
//			InputIterator last) :
//			T_(T), Q_(first, last) {
//	}

	template<typename KnotInputIterator, typename VectorInputIterator>
	bspline(KnotInputIterator T_first, KnotInputIterator T_last,
			VectorInputIterator Q_first, VectorInputIterator Q_last) :
			T_(T_first, T_last), Q_(Q_first, Q_last) {

	}

	~bspline() {

	}

	/**
	 * @brief Returns the degree of this B-spline.
	 * @return
	 */
	size_t degree() const {
		return this->compute_degree_(this->T_.size(), this->Q_.size());
	}

	std::pair<parameter, parameter> domain() const {
		return std::make_pair(this->T_[this->degree()],
				this->T_[this->Q_.size()]);
	}

	/**
	 * @brief Returns the knot vector.
	 * @return
	 */
	const knotvector_type& knot_vector() const {
		return this->T_;
	}

	/**
	 * @brief Returns the control points.
	 * @return
	 */
	const control_points& controls() const {
		return this->Q_;
	}

	/**
	 * @brief Returns the mutable control points.
	 * @return
	 */
	control_points& controls() {
		return this->Q_;
	}

	void insert(const parameter& t) {
		this->insert_knot_(t);
	}

	/**
	 * @brief Subdivides this at a parameter @a t.
	 *
	 * @param t Parameter to subdivide.
	 * @param selection
	 *
	 * @return Another.
	 */
	bspline subdivide(const parameter& t, gkselection selection = GK::Upper) {
		const gksize degree = this->compute_degree_(this->T_.size(),
				this->Q_.size());
		const gksize order = degree + 1;

		if (t < this->T_[degree] || t > this->T_[this->Q_.size()]) {
			return bspline();
		}

		for (gksize i = 0; i < order; ++i) {
			this->insert_knot_(t);
		}

		typename knotvector_type::const_iterator first = this->T_.begin()
				+ order;
		typename knotvector_type::const_iterator end = this->T_.begin()
				+ this->Q_.size() /*+ order*/;

		typename knotvector<parameter>::const_iterator position =
				std::upper_bound(first, end, t);

		const gksize k = std::distance(this->T_.begin(), position);

		if (selection == GK::Upper) {
			const knotvector<parameter> other_knotvector(position - order,
					this->T_.end());
			const control_points other_control_points(
					this->Q_.begin() + k - order, this->Q_.end());

			this->T_.erase(position, this->T_.end());
			this->Q_.erase(this->Q_.begin() + k - order, this->Q_.end());

			return bspline(other_knotvector, other_control_points);

		} else {
			const knotvector<parameter> other_knotvector(this->T_.begin(),
					position);
			const control_points other_control_points(this->Q_.begin(),
					this->Q_.begin() + k - order);

			this->T_.erase(this->T_.begin(), position - order);
			this->Q_.erase(this->Q_.begin(), this->Q_.begin() + k - order);

			return bspline(other_knotvector, other_control_points);
		}
	}

	/**
	 * @brief
	 * @param selection
	 * @return
	 */
	bspline subdivide(gkselection selection = GK::Upper) {
		const parameter t_max = this->T_[this->Q_.size()];
		const parameter t_min = this->T_[this->compute_degree_(this->T_.size(),
				this->Q_.size())];

		return subdivide(0.5 * (t_max + t_min), selection);
	}

//	void derivatise() {
//		const size_t degree = this->compute_degree_(this->T_.size(),
//				this->Q_.size());
//		const parameter Zero = parameter(GK_FLOAT_ZERO);
//
//		for (size_t i = 0; i < this->Q_.size() - 1; ++i) {
//			const parameter denominator = this->T_[i + degree + 1]
//					- this->T_[i + 1];
//			this->Q_[i] =
//					(denominator == Zero) ?
//							Zero :
//							degree * (this->Q_[i + 1] - this->Q_[i])
//									/ denominator;
//		}
//
//		this->Q_.pop_back();
//
//		this->T_.pop_back();
//		this->T_.erase(this->T_.begin());
//	}

	vector_type operator()(const parameter& t) const {
		std::vector<parameter> N(this->Q_.size());
		basis_function(compute_degree_(this->T_.size(), this->Q_.size()),
				this->T_.begin(), this->T_.end(), t, N.begin());

//		vector_type r; //= N[0] * this->Q_[0];
//		for (size_t j = 0; j < this->Q_.size(); j++) {
//			r += N[j] * this->Q_[j];
//		}

		return std::inner_product(this->Q_.begin(), this->Q_.end(), N.begin(),
				vector_type());
//		return r;
	}

	template<typename InputIterator, typename OutputIterator>
	OutputIterator operator()(InputIterator first, InputIterator last,
			OutputIterator result) const {

		return result;
	}

	bspline& operator=(const bspline& rhs) {
		if (&rhs == this) {
			return *this;
		}

		this->T_ = rhs.T_;
		this->Q_ = rhs.Q_;

		return *this;
	}

private:
	knotvector_type T_;
	control_points Q_;

private:
	void insert_knot_(const parameter& t) {
		const size_t degree = this->compute_degree_(this->T_.size(),
				this->Q_.size());

//		const typename knotvector<parameter>::const_iterator valid_first =
//				this->T_.begin() + degree;
//		const typename knotvector<parameter>::const_iterator valid_last =
//				this->T_.begin() + this->Q_.size();

		if (t < this->T_[degree] || t > this->T_[this->Q_.size()]) {
			return;
		}

		// section number
		const size_t k = knot_segment(this->T_, degree, t);

		this->Q_.insert(this->Q_.begin() + k, this->Q_[k]);

		for (size_t j = k; j > k - degree; --j) {
			const parameter denominator = this->T_[j + degree] - this->T_[j];
			const parameter alpha =
					(denominator == parameter(GK_FLOAT_ZERO)) ?
							parameter(GK_FLOAT_ZERO) :
							(t - this->T_[j]) / (denominator);

			this->Q_[j] = alpha * this->Q_[j]
					+ (parameter(GK_FLOAT_ONE) - alpha) * this->Q_[j - 1];
		}

		this->T_.insert(t);
	}

};

/**
 * @brief Computes the parameter domain of the B-spline @a r.
 * @param r
 * @return
 */
template<typename Vector>
std::pair<typename bspline<Vector>::parameter,
		typename bspline<Vector>::parameter> domain(const bspline<Vector>& r) {
	const typename bspline<Vector>::knotvector_type& T = r.knot_vector();
	return std::make_pair(T[r.degree()], T[r.controls().size()]);
}

template<typename Vector>
std::pair<bspline<Vector>, bspline<Vector> > subdivide(const bspline<Vector>& x,
		const typename bspline<Vector>::parameter& t) {
	bspline<Vector> y = x;
	const bspline<Vector> z = y.subdivide(t, GK::Upper);
	return std::make_pair(y, z);
}

template<typename Vector>
bspline<Vector> subdivide(const bspline<Vector>& x,
		const typename bspline<Vector>::parameter& a,
		const typename bspline<Vector>::parameter& b) {
	bspline<Vector> y = x;
	y.subdivide(a, GK::Lower);
	y.subdivide(b, GK::Upper);

	return y;
}

template<typename Vector, typename Parameter>
segment<Vector> linear(const bspline<Vector>& x,
		typename vector_traits<Vector>::value_type& max_distance) {

	max_distance = typename vector_traits<Vector>::value_type(GK_FLOAT_ZERO);
	typedef typename bspline<Vector>::parameter parameter;
	std::pair<parameter, parameter> D = x.domain();
	bspline<Vector> y = x;
	y.subdivide(D.first, GK::Lower);
	y.subdivide(D.second, GK::Upper);

	const segment<Vector> l(y.controls().front(), y.controls().back());
	typedef typename bspline<Vector>::control_points controls;
	for (typename controls::iterator p = y.controls().begin();
			p != y.controls().end(); ++p) {
		const typename geometry_traits<segment<Vector> >::parameter t = nearest(
				l, *p);

		const typename vector_traits<Vector>::value_type d = norm(l(t) - *p);
		max_distance = std::max(max_distance, d);
	}
	return l;
}

template<typename Vector>
aabb<Vector> boundary(const bspline<Vector>& x) {
	return aabb<Vector>(x.controls().begin(), x.controls().end());
}

template<typename Vector>
struct curve_traits<bspline<Vector> > : public geometry_traits<bspline<Vector> > {
	typedef geometry_traits<bspline<Vector> > base;
	typedef typename base::parameter parameter;
	typedef typename vector_traits<Vector>::value_type value_type;
	typedef value_type distance_type;
	typedef value_type length_type;
	typedef aabb<Vector> boundary_type;

	static std::pair<parameter, parameter> domain(const bspline<Vector>& x) {
		return x.domain();
	}

	static boundary_type boundary(const bspline<Vector>& x) {
		return boundary_type(x.controls().begin(), x.controls().end());
	}

};

template<typename Vector, typename Parameter>
bspline<Vector, Parameter> derivatatise(const bspline<Vector, Parameter>& r) {

	const Parameter Zero(GK_FLOAT_ZERO);

//	typedef typename division_result<Vector, Parameter>::value_type dV_t;

	const size_t degree = r.degree();
	typedef typename bspline<Vector, Parameter>::knotvector_type knotvector;
	const knotvector T = r.knot_vector();
	const typename bspline<Vector, Parameter>::control_points Q = r.controls();

	typename bspline<Vector, Parameter>::control_points P(Q.size() - 1);
	for (size_t i = 0; i < P.size(); ++i) {
		const Parameter dt = T[i + degree + 1] - T[i + 1];
		P[i] = (dt == Zero) ? Zero : degree * (Q[i + 1] - Q[i]) / dt;
	}

	typename knotvector::iterator U_first = T.begin();
	std::advance(U_first, 1);

	typename knotvector::iterator U_last = T.end();
	std::advance(U_last, -1);

	return bspline<Vector, Parameter>(U_first, U_last, P.begin(), P.end());
}

}  // namespace gk

///**
// * @brief B-Spline (Basis spline).
// *
// * @author Takuya Makimoto
// * @date 2015
// */
//template<GK::DimensionNumber Size, typename Parameter, typename Scalar>
//class gkbspline {
//public:
//	typedef vectors<Size, Scalar> controls_type;
//	typedef Parameter parameter;
//	typedef vector<Size, Scalar> dependence;
////	typedef typename controls_type::iterator iterator;
////	typedef typename controls_type::const_iterator const_iterator;
////	typedef typename controls_type::reverse_iterator reverse_iterator;
////	typedef typename controls_type::const_reverse_iterator const_reverse_iterator;
//
//	typedef enum {
//		InsertControl, InsertOrder,
//	} InsertType;
//
//	static const parameter Start;
//	static const parameter End;
//
//private:
//	static const gkuint InvalidOrder_;
//
//public:
//	/**
//	 * @brief Default constructor.
//	 * The object created by this method is invalid.
//	 */
//	gkbspline() :
//			controls_(), T_() {
//
//	}
//
//	/**
//	 * @brief Copy constructor.
//	 */
//	gkbspline(const gkbspline& other) :
//			controls_(other.controls_), T_(other.T_) {
//
//	}
//
//	gkbspline(const controls_type& controls, const gkknotvector& T) :
//			controls_(controls), T_(T) {
//
//	}
//
//	/**
//	 * @brief Destructor.
//	 */
//	~gkbspline() {
//
//	}
//
////	/**
////	 * @copydoc gkvector::valid()
////	 */
////	bool valid() const;
//
//	gksize order() const;
//
//	const controls_type& controls() const {
//		return this->controls_;
//	}
//
//	/**
//	 * @brief
//	 */
//	void set_control_point(gkuint index, const gkvector<Size>& position) {
//		this->controls_[index] = position;
//	}
//
//	/**
//	 * @brief Returns the constant reference of the knot vector object.
//	 */
//	const gkknotvector& knotvector() const {
//		return this->T_;
//	}
//
////	template<class Generator>
////	void distribute_knot(Generator generator) {
////		this->T_.distribute(generator);
//////        normalize_knotvector_();
////	}
//
//	/**
//	 * @brief Inserts a knot.
//	 * @param t Knot to insert.
//	 * @param value_type
//	 * @return
//	 */
//	void insert_knot(parameter t, InsertType value_type = gkbspline::InsertControl);
//
//	/**
//	 *
//	 * @return
//	 */
//	std::pair<gkknotvector::const_iterator, gkknotvector::const_iterator> valid_range() const;
//
//	/**
//	 * @brief Subdivides at @a t.
//	 * @param t Parameter at subdividing position.
//	 * @return
//	 */
//	void subdivide(parameter t, gkselection selection = GK::Upper);
//
//	/**
//	 *
//	 */
//	void differentiate();
//
//	/**
//	 * @brief Elevates the order.
//	 */
//	void elevate();
//
//	/**
//	 * @brief Reduces the order.
//	 */
//	void reduce();
//
////	void translate(const gkvector<DimensionNumber>& translation);
//
////	void rotate(const gkrotation<DimensionNumber>& rotation);
//
////	void scale(gkfloat factor);
//
//	/**
//	 * @brief Computes the position at the parameter @a t.
//	 * @param t
//	 * @return
//	 */
//	vector<Size, Scalar> operator()(const parameter& t) const {
//
//	}
//
////	const gkvector<DimensionNumber>& operator[](gkuint index) const;
////	gkvector<DimensionNumber>& operator[](gkuint index);
//
//	/**
//	 *
//	 */
//	gkbspline& operator=(const gkbspline& rhs) {
//		if (&rhs == this) {
//			return *this;
//		}
//
//		this->controls_ = rhs.controls_;
//		this->T_ = rhs.T_;
//
//		return *this;
//	}
//
//private:
//	controls_type controls_;
//	gkknotvector T_;
//
//private:
//	/**
//	 * @brief Computes the order of this spline.
//	 * @return
//	 */
//	gkuint order_() const;
//
//	void normalize_knotvector_();
//
//	void insert_knot_inserted_control_(parameter t);
//
//	void insert_knots_(const gkknotvector& U);
//};
//
//}  // namespace gk

//template<GK::DimensionNumber DimensionNumber>
//std::pair<gkfloat, gkfloat> parameter_range(
//		const gkbspline<DimensionNumber>& bspline) {
//	const gkknotvector& T = bspline.knotvector();
//
//	return std::make_pair(T[bspline.order() - 1], T[bspline.controls().size()]);
//}
//
///**
// * @todo Implementation.
// * @param bspline
// * @return
// */
//template<GK::DimensionNumber DimensionNumber>
//gklength length(const gkbspline<DimensionNumber>& bspline) {
//	return 0.0;
//}
//
//template<GK::DimensionNumber DimensionNumber>
//gkaabb<DimensionNumber> boundary(const gkbspline<DimensionNumber>& bspline) {
//	const std::pair<gkfloat, gkfloat> t = parameter_range(bspline);
//	gkbspline<DimensionNumber> r = bspline;
//
//	r.subdivide(t.first, GK::Lower);
//	r.subdivide(t.second, GK::Upper);
//
//	return boundary(r.controls());
//}
//
//template<GK::DimensionNumber DimensionNumber>
//std::pair<gkbspline<DimensionNumber>, gkbspline<DimensionNumber> > subdivide(
//		const gkbspline<DimensionNumber>& bspline, gkfloat t) {
//	gkbspline<DimensionNumber> upper = bspline;
//	gkbspline<DimensionNumber> lower = bspline;
//
//	upper.subdivide(t, GK::Upper);
//	lower.subdivide(t, GK::Lower);
//
//	return std::make_pair(upper, lower);
//}
//
//template<GK::DimensionNumber DimensionNumber>
//class const_gkbspline {
//public:
//	static const gkbspline<DimensionNumber>& Invalid();
//};
//
//template<GK::DimensionNumber DimensionNumber, template<DimensionNumber> class Curve>
//gkbspline<DimensionNumber> convert_bspline(const Curve<DimensionNumber>& curve);

#endif /* BSPLINE_BSPLINE_H_ */
