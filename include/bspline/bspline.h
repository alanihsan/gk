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
#include "basis.h"

namespace gk {

template<typename KnotVector>
struct knotvector_traits {
	typedef typename KnotVector::value_type value_type;
	typedef typename KnotVector::const_reference const_reference;
	typedef typename KnotVector::const_iterator const_iterator;
	typedef typename KnotVector::const_reverse_iterator const_reverse_iterator;
};

template<typename KnotVector>
typename knotvector_traits<KnotVector>::const_iterator begin(
		const KnotVector& T);
template<typename KnotVector>
typename knotvector_traits<KnotVector>::const_iterator end(const KnotVector& T);
template<typename KnotVector>
typename knotvector_traits<KnotVector>::const_reverse_iterator rbegin(
		const KnotVector& T);
template<typename KnotVector>
typename knotvector_traits<KnotVector>::const_reverse_iterator rend(
		const KnotVector& T);

template<typename KnotVector>
size_t size_of(const KnotVector& T);

template<typename KnotVector>
std::pair<typename knotvector_traits<KnotVector>::value_type,
		typename knotvector_traits<KnotVector>::value_type> domain(
		size_t degree, const KnotVector& T) {

	return std::make_pair(T[degree], T[size_of(T) - degree + 1]);
}

template<typename KnotVector>
typename knotvector_traits<KnotVector>::const_iterator erase_element(
		KnotVector& T,
		typename knotvector_traits<KnotVector>::const_iterator position) {
	return T.erase(position);
}

template<typename KnotVector>
typename knotvector_traits<KnotVector>::const_iterator erase_elements(
		KnotVector& T,
		typename knotvector_traits<KnotVector>::const_iterator first,
		typename knotvector_traits<KnotVector>::const_iterator last) {
	return T.erase(first, last);
}

/**
 * @brief B-spline (Basis spline).
 *
 * @tparam Vector Type of a control point.
 *
 * @author Takuya Makimoto
 * @date 2015
 */
template<typename Vector, typename Parameter>
class bspline: public curve<free_curve_tag, Vector> {
public:
	typedef Vector vector_type;
	typedef std::vector<Parameter> knotvector_type;
	typedef std::vector<Vector> control_points;

private:
	/**
	 * @brief Computes the degree of this B-spline.
	 * @param knotvector_size
	 * @param controls_size
	 * @return
	 */
	static size_t compute_degree_(size_t knotvector_size,
			size_t controls_size) {
		return knotvector_size - controls_size - 1;
	}

public:

	/**
	 * @brief Default constructor.
	 */
	bspline() :
			T_(), Q_() {
	}

	/**
	 * @brief Copy constructor.
	 * @param other
	 */
	bspline(const bspline& other) :
			T_(other.T_), Q_(other.Q_) {

	}

	/**
	 *
	 * @param T A knot vector.
	 * @param Q An object of control points.
	 */
	bspline(const knotvector_type& T, const control_points& Q) :
			T_(T), Q_(Q) {
	}

	/**
	 * @brief
	 * @tparam KnotInputIterator
	 * @tparam VectorInputIterator
	 * @param T_first
	 * @param T_last
	 * @param Q_first
	 * @param Q_last
	 */
	template<typename KnotInputIterator, typename VectorInputIterator>
	bspline(KnotInputIterator T_first, KnotInputIterator T_last,
			VectorInputIterator Q_first, VectorInputIterator Q_last) :
			T_(T_first, T_last), Q_(Q_first, Q_last) {

	}

	/**
	 * @brief Destructor.
	 */
	~bspline() {
	}

	/**
	 * @brief Returns the degree of this B-spline.
	 * @return
	 */
	size_t degree() const {
		return this->compute_degree_(size_of(this->T_), this->Q_.size());
	}

//	std::pair<parameter, parameter> domain() const {
//		return std::make_pair(this->T_[this->degree()],
//				this->T_[this->Q_.size()]);
//	}

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

	template<typename Parameter>
	void insert(const Parameter& t) {
		this->insert_knot_(t);
	}

	/**
	 * @brief Subdivides this at a parameter @a t.
	 *
	 * @param t Parameter to subdivide.
	 * @param selection
	 *
	 * @return The other.
	 */
	template<typename Parameter>
	bspline subdivide(const Parameter& t, gkselection selection = GK::Upper) {
		const gksize degree = this->compute_degree_(this->T_.size(),
				this->Q_.size());
		const gksize order = degree + 1;

		if (t < this->T_[degree] || t > this->T_[this->Q_.size()]) {
			return bspline();
		}

		for (gksize i = 0; i < order; ++i) {
			this->insert_knot_(t);
		}

		typedef typename knotvector_traits<knotvector_type>::const_iterator T_const_iterator;
		T_const_iterator first = this->T_.begin() + order;
		T_const_iterator end = this->T_.begin() + this->Q_.size() /*+ order*/;

		T_const_iterator position = std::upper_bound(first, end, t);

		const gksize k = std::distance(this->T_.begin(), position);

		if (selection == GK::Upper) {
			const knotvector_type other_knotvector(position - order,
					this->T_.end());
			const control_points other_control_points(
					this->Q_.begin() + k - order, this->Q_.end());

			erase_elements(this->T_, position, this->T_.end());
			this->Q_.erase(this->Q_.begin() + k - order, this->Q_.end());

			return bspline(other_knotvector, other_control_points);

		} else {
			const knotvector_type other_knotvector(this->T_.begin(), position);
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
		typedef typename knotvector_traits<knotvector_type>::value_type parameter;
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

	Vector operator()(const Parameter& t) const {
		std::vector<Parameter> N(this->Q_.size());
		basis_function(compute_degree_(this->T_.size(), this->Q_.size()),
				this->T_.begin(), this->T_.end(), t, N.begin());

		return std::inner_product(this->Q_.begin(), this->Q_.end(), N.begin(),
				Vector());
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
	template<typename Parameter>
	void insert_knot_(const Parameter& t) {
		const size_t degree = this->compute_degree_(this->T_.size(),
				this->Q_.size());

//		const typename knotvector<parameter>::const_iterator valid_first =
//				this->T_.begin() + degree;
//		const typename knotvector<parameter>::const_iterator valid_last =
//				this->T_.begin() + this->Q_.size();

		const std::pair<Parameter, Parameter> D = domain(degree, this->T_);
		if (t < D.first || t > D.second) {
			return;
		}

		// section number
		const size_t k = knot_segment(this->T_, degree, t);

		this->Q_.insert(this->Q_.begin() + k, this->Q_[k]);

		for (size_t j = k; j > k - degree; --j) {
			const Parameter denominator = this->T_[j + degree] - this->T_[j];
			const Parameter alpha =
					(denominator == Parameter(GK_FLOAT_ZERO)) ?
							Parameter(GK_FLOAT_ZERO) :
							(t - this->T_[j]) / (denominator);

			this->Q_[j] = alpha * this->Q_[j]
					+ (Parameter(GK_FLOAT_ONE) - alpha) * this->Q_[j - 1];
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
	return domain(r.degree(), r.knot_vector());
}

template<typename Vector, typename Parameter>
std::pair<bspline<Vector>, bspline<Vector> > subdivide(const bspline<Vector>& r,
		const Parameter& t) {
	bspline<Vector> p = r;
	const bspline<Vector> q = p.subdivide(t, GK::Upper);
	return std::make_pair(p, q);
}

template<typename Vector, typename Parameter>
bspline<Vector> subdivide(const bspline<Vector>& x, const Parameter& a,
		const Parameter& b) {
	bspline<Vector> y = x;
	y.subdivide(a, GK::Lower);
	y.subdivide(b, GK::Upper);

	return y;
}

template<typename Vector, typename KnotVector, typename Parameter,
		typename OutputIterator>
OutputIterator subdivide(const bspline<Vector, KnotVector>& x,
		const Parameter& t, OutputIterator result) {
	bspline<Vector> y = x;
	*result = y.subdivide(t, GK::Lower);
	++result;
	*result = y;
	return result;
}

template<typename Vector, typename KnotVector>
std::pair<Vector, Vector> linearize(const bspline<Vector, KnotVector>& x,
		typename vector_traits<Vector>::value_type& max_distance) {

	max_distance = typename vector_traits<Vector>::value_type(
	GK_FLOAT_ZERO);
	typedef typename bspline<Vector>::parameter parameter;
	std::pair<parameter, parameter> D = x.domain();
	bspline<Vector> y = x;
	y.subdivide(D.first, GK::Lower);
	y.subdivide(D.second, GK::Upper);

	std::pair<Vector, Vector> result = std::make_pair(y.controls().front(),
			y.controls().back());
	typedef typename bspline<Vector>::control_points controls;
	for (typename controls::iterator p = y.controls().begin();
			p != y.controls().end(); ++p) {
		const Vector v = alg::nearest_to_line(result.first,
				normalize(result.second - result.first), *p);

		const typename vector_traits<Vector>::value_type d = norm(v - *p);
		max_distance = std::max(max_distance, d);
	}

	return result;
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

template<typename Vector>
bspline<Vector> derivatatise(const bspline<Vector>& r) {

	typedef typename bspline<Vector>::knotvector_type knotvector;
	typedef typename knotvector_traits<knotvector>::value_type parameter;
	const parameter Zero(GK_FLOAT_ZERO);

//	typedef typename division_result<Vector, parameter>::value_type dV_t;

	const size_t degree = r.degree();
	const knotvector T = r.knot_vector();
	const typename bspline<Vector, parameter>::control_points Q = r.controls();

	typename bspline<Vector, parameter>::control_points P(Q.size() - 1);
	for (size_t i = 0; i < P.size(); ++i) {
		const parameter dt = T[i + degree + 1] - T[i + 1];
		P[i] = (dt == Zero) ? Zero : degree * (Q[i + 1] - Q[i]) / dt;
	}

	typedef typename knotvector_traits<knotvector>::const_iterator T_const_iterator;
	T_const_iterator U_first = T.begin();
	std::advance(U_first, 1);

	T_const_iterator U_last = T.end();
	std::advance(U_last, -1);

	return bspline<Vector>(U_first, U_last, P.begin(), P.end());
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
