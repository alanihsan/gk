/*
 * line.h
 *
 *  Created on: 2014/07/16
 *      Author: makitaku
 */

#ifndef PRIMITIVE_LINE_H_
#define PRIMITIVE_LINE_H_

#include "../gkvector.h"
#include "../gkcurve.h"
#include "../gkaabb.h"

#include <vector>

#include <Eigen/Core>

namespace gk {

template<std::size_t Dimension, typename T>
struct vector_t {
	typedef void type;
	template<typename Vector>
	static type assign(const Vector& v);
};

template<typename T>
struct vector_t<GK::GK_2D, T> {
	typedef Eigen::Matrix<T, 2, GK::GK_2D> type;

	template<typename Vector>
	static type assign(const Vector& v) {
		type r;
		r[GK::X] = v[GK::X];
		r[GK::Y] = v[GK::Y];
		return r;
	}
};

template<typename T>
struct vector_t<GK::GK_3D, T> {
	typedef Eigen::Matrix<T, 1, GK::GK_3D> type;

	template<typename Vector>
	static type assign(const Vector& v) {
		type r;
		r[GK::X] = v[GK::X];
		r[GK::Y] = v[GK::Y];
		r[GK::Z] = v[GK::Z];

		return r;
	}
};

/**
 * @brief Line.
 * This object computes with the geometric formula of lines as below.@n
 * @f{displaymath}{
 * 	\mathbf{r}(t) = \mathbf{p} + t\mathbf{d}
 * @f}@n
 * where @f$\mathbf{r}@f$ is line function, @f$\mathbf{p}@f$ is the reference point,
 * @f$\mathbf{d}@f$ is the direction (unit vector) and @f$t@f$ is the parameter.
 * @n
 * This object has a reference point @f$\mathbf{p}@f$ and a direction
 * @f$\mathbf{d}@f$ which are inner vector.
 *
 * @author Takuya Makimoto
 * @date 2015/12/01
 */
template<std::size_t DimensionSize, typename T>
class xline/*: public geometry<line_tag, Vector>*/{
public:
//	typedef Vector vector_type;
	static const std::size_t Dimension = DimensionSize;
	typedef direction<Dimension> direction_type;

private:
	typedef typename vector_t<DimensionSize, T>::type vector_type;

public:
	/**
	 * @brief Default contsruction.
	 */
	xline() :
			ref_(), direction_() {
	}

	/**
	 * @brief Copy constructor.
	 * @param other An other object.
	 */
	xline(const xline& other) :
			ref_(other.ref_), direction_(other) {

	}

	template<typename Position>
	xline(const Position& reference, const direction_type& direction) :
			ref_(vector_t<DimensionSize, T>::assign(reference)), direction_(
					direction) {
	}

	template<typename Position>
	xline(const Position& start, const Position& end) :
			ref_(vector_t<DimensionSize, T>::assign(start)), direction_(
					vector_traits<vector_type>::begin(end - start)) {
	}

//	line(const std::pair<vector_type, vector_type>& pair) :
//			ref_(pair.first), direction_(pair.second - pair.first) {
//	}

	~xline() {
	}

	point<DimensionSize, T> reference() const {
		return this->ref_;
	}

//	vector_type& reference() {
//		return this->ref_;
//	}
	template<typename Position>
	void reference(const Position& r) {
		for (std::size_t i = 0; i < DimensionSize; ++i) {
			this->ref_[i] = r[i];
		}
	}

	/**
	 *
	 * @return
	 */
	const direction_type& tangent_direction() const {
		return this->direction_;
	}

	/**
	 *
	 * @param d
	 */
	void tangent_direction(const direction_type& d) {
		this->direction_ = d;
	}

	/**
	 * @brief Computes the position at the parameter @a t.
	 *
	 * @f[
	 * 		\mathbf{r}(t) = \mathbf{p} + t\mathbf{d}
	 * @f]
	 *
	 * @tparam Parameter Type of a parameter.
	 * @param t Parameter.
	 * @return Position vector.
	 */
	template<typename Parameter>
	point<DimensionSize, T> operator()(const Parameter& t) const {
		return this->ref_ + t * this->direction_;
	}

	xline& operator=(const xline& m) {
		if (&m == this) {
			return *this;
		}

		this->ref_ = m.ref_;
		this->direction_ = m.direction_;

		return *this;
	}

private:
	vector_type ref_;
	direction_type direction_;
};

template<std::size_t DimensionSize, typename T>
struct geometry_traits<xline<DimensionSize, T> > {
	typedef line_tag category;
	static const std::size_t Dimension = DimensionSize;
};

//template<typename Vector, bool Upper = true>
//class ray: public curve<Vector, ray_tag,
//		typename vector_traits<Vector>::value_type> {
//public:
//	ray() :
//			ref_(), delta_() {
//	}
//
//	ray(const ray& other) :
//			ref_(other.ref_), delta_(other.delta_) {
//	}
//
//	~ray() {
//	}
//
//	vector_type operator()(const parameter& t) const {
//		const typename norm_result<vector_type>::value_type L = norm(this->delta_);
//		return this->ref_ + t * this->delta_ / L;
//	}
//
//private:
//	vector_type ref_;
//	vector_type delta_;
//};
//
/**
 * @brief Segment.
 *  @f[
 * \mathbf{l}(t) = \mathbf{l}(\mathbf{m}(s + t \times (e - s)))
 * @f]
 * @f$\mathbf{l}@f$ : Segment.@n
 * @f$t@f$ : Parameter.@n
 * @f$\mathbf{m}@f$ : Line.@n
 * @f$s@f$ : Start position parameter on the line, @a m.@n
 * @f$e@f$ : End position parameter on the line, @a m.@n
 *
 * @author Takuya Makimoto
 * @date 2015/12/01
 */
template<std::size_t DimensionSize,typename T>
class segment{
public:
//	typedef typename vector_traits<Vector>::value_type value_type;


private:
	typedef vector_t<DimensionSize,T>::type vector_type;
public:
	segment() :
			edge_() {
	}

	segment(const segment& other) :
			edge_() {
		this->edge_[GK::StartEdge] = other.edge_[GK::StartEdge];
		this->edge_[GK::EndEdge] = other.edge_[GK::EndEdge];
	}

	segment(const vector_type& start, const vector_type& end) :
			edge_() {
		this->edge_[GK::StartEdge] = start;
		this->edge_[GK::EndEdge] = end;
	}

	segment(const std::pair<vector_type, vector_type>& pair) :
			edge_() {
		this->edge_[GK::StartEdge] = pair.first;
		this->edge_[GK::EndEdge] = pair.second;
	}

	~segment() {
	}

	void inverse() {
		std::swap(this->edge_[GK::StartEdge], this->edge_[GK::EndEdge]);
	}

	const vector_type& start() const {
		return this->edge_[GK::StartEdge];
	}

	void start(const vector_type& start) {
		this->edge_[GK::StartEdge] = start;
	}

	const vector_type& end() const {
		return this->edge_[GK::EndEdge];
	}

	void end(const vector_type& end) {
		this->edge_[GK::EndEdge] = end;
	}

	const vector_type& operator[](size_t index) const {
		return this->edge_[index];
	}

	vector_type& operator[](size_t index) {
		return this->edge_[index];
	}

	template<typename Parameter>
	point<DimensionSize,T> operator()(const Parameter& t) const {
		const vector_type v = this->edge_[GK::EndEdge]
				- this->edge_[GK::StartEdge];
		return this->edge_[GK::StartEdge] + t * v;
	}

	segment& operator=(const segment& rhs) {
		if (&rhs == this) {
			return *this;
		}

		this->edge_[GK::StartEdge] = rhs.edge_[GK::StartEdge];
		this->edge_[GK::EndEdge] = rhs.edge_[GK::EndEdge];

		return *this;
	}

private:
	vector_type edge_[GK::EdgeSize];
};

//template<typename Vector>
//struct curve_traits<segment<Vector> > {
//	typedef typename geometry_traits<segment<Vector> >::parameter parameter;
//	typedef aabb<Vector> boundary_type;
//
//	static std::pair<parameter, parameter> domain(const segment<Vector>&) {
//		return std::make_pair(parameter(GK_FLOAT_ZERO), parameter(GK_FLOAT_ONE));
//	}
//
//};

template<typename Vector>
direction<vector_traits<Vector>::Dimension> direction_of(
		const line<Vector>& l) {
	return l.tangent_direction();
}

template<typename Vector>
direction<vector_traits<Vector>::Dimension> direction_of(
		const segment<Vector>& l) {
	return direction_of(l[GK::EndEdge] - l[GK::StartEdge]);
}

template<typename Vector>
typename vector_traits<Vector>::value_type length(const segment<Vector>& x) {
	return norm(x.start() - x.end());
}

template<typename Vector>
std::pair<typename segment<Vector>::parameter,
		typename segment<Vector>::parameter> domain(const segment<Vector>&) {
	typedef typename segment<Vector>::parameter parameter;
	return std::make_pair(parameter(GK_FLOAT_ZERO), parameter(GK_FLOAT_ONE));
}

template<typename Vector>
aabb<Vector> boundary(const segment<Vector>& x) {
	return aabb<Vector>(x.start(), x.end());
}

template<typename Vector>
class polyline/*: public curve<curve_tag, Vector>*/{
public:
	typedef std::vector<Vector> container_type;
	typedef Vector vector_type;
	typedef typename vector_traits<Vector>::value_type value_type;
//	typedef Parameter parameter;

	typedef typename container_type::iterator iterator;
	typedef typename container_type::const_iterator const_iterator;
	typedef typename container_type::reverse_iterator reverse_iterator;
	typedef typename container_type::const_reverse_iterator const_reverse_iterator;

private:
	template<typename InputIterator, typename OutputIterator>
	static OutputIterator lengths(InputIterator first, InputIterator last,
			OutputIterator result) {

		*result = vector_traits<Vector>::value_type(GK_FLOAT_ZERO);
		++result;
		while (first != last) {
			*result = vector_traits<Vector>::norm(*first++ - *first);
			++result;
		}

		return result;
	}

public:
	polyline() :
			X_() {

	}

	polyline(const polyline& other) :
			X_(other.X_) {

	}

	template<typename InputIterator>
	polyline(InputIterator first, InputIterator last) :
			X_(first, last) {

	}

	~polyline() {

	}

	bool empty() const {
		return this->X_.empty();
	}

	const_iterator begin() const {
		return this->X_.begin();
	}

	iterator begin() {
		return this->X_.begin();
	}

	const_iterator end() const {
		return this->X_.end();
	}

	iterator end() {
		return this->X_.end();
	}

	const_reverse_iterator rbegin() const {
		return this->X_.rbegin();
	}

	reverse_iterator rbegin() {
		return this->X_.rbegin();
	}

	const_reverse_iterator rend() const {
		return this->X_.rend();
	}

	reverse_iterator rend() {
		return this->X_.rend();
	}

	template<typename Parameter>
	polyline subdivide(const Parameter& t, gkselection select = GK::Upper) {

	}

	const vector_type& operator[](size_t index) const {
		return this->X_[index];
	}

	vector_type& operator[](size_t index) {
		return this->X_[index];
	}

	template<typename Parameter>
	vector_type operator()(const Parameter& t) const {
		std::vector<value_type> L;
		L.reserve(this->X_.size());

		lengths(this->X_.begin(), this->X_.end(), std::inserter(L, L.begin()));

		const value_type p = t * L.back();
		typename std::vector<value_type>::iterator bound = std::upper_bound(
				L.begin(), L.end(), p);

		const size_t n = std::distance(L.begin(), bound);
		const Parameter u = (*bound - p) / L.back();

		return (Parameter(GK_FLOAT_ONE) - u) * this->X_[n - 1] + u * this->X_[n];
	}

	polyline& operator=(const polyline& rhs) {
		if (&rhs == this) {
			return *this;
		}

		this->X_ = rhs.X_;
		return *this;
	}

private:
	container_type X_;
};

} // namespace gk

#endif /* PRIMITIVE_LINE_H_ */
