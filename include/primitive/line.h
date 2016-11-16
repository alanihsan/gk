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

namespace gk {

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
template<typename T, std::size_t Dimension = GK::GK_3D>
<<<<<<< HEAD
class line: public geometry<line_tag, T, Dimension> {
=======
class line: public geometry<line_tag> {
>>>>>>> 5ed5cf2db97b02bd867d82a1513c790b646cb925
public:
//	typedef T value_type;

	typedef typename vector_type<T, Dimension>::type vector_type;
	typedef direction<Dimension> direction_type;

private:
	template<typename Vector>
	static vector_type Reference_(const Vector& v, direction_tag<GK::GK_2D>) {
		return vector_type(v[GK::X], v[GK::Y]);
	}

	template<typename Vector>
	static vector_type Reference_(const Vector& v, direction_tag<GK::GK_3D>) {
		return vector_type(v[GK::X], v[GK::Y], v[GK::Z]);
	}

public:
	/**
	 * @brief Default contsruction.
	 */
	line() :
			ref_(), direction_() {
	}

	/**
	 * @brief Copy constructor.
	 * @param other An other object.
	 */
	line(const line& other) :
			ref_(other.ref_), direction_(other) {

	}

	line(const vector_type& reference, const direction_type& direction) :
			ref_(reference), direction_(direction) {
	}

	line(const vector_type& start, const vector_type& end) :
			ref_(start), direction_(
					vector_traits<vector_type>::begin(end - start)) {
	}

	template<typename Vector>
	line(const Vector& reference, const direction_type& direction) :
			ref_(line::Reference_(reference, dimension_tag<Dimension>())), direction_(
					direction) {
	}

	template<typename Vector>
	line(const Vector& start, const Vector& end) :
			ref_(line::Reference_(start, dimension_tag<Dimension>())), direction_() {
	}

	~line() {
	}

	const vector_type& reference() const {
		return this->ref_;
	}

	void reference(const vector_type& r) {
		this->ref_ = r;
	}

	template<typename Vector>
	void reference(const Vector& r) {
		this->ref_ = r;
	}

	/**
	 * @brief Returns the direction.
	 * @return
	 */
	const direction_type& direction() const {
		return this->direction_;
	}

	/**
	 * @brief Sets a direction.
	 * @param d
	 */
	void direction(const direction_type& d) {
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
	vector_type operator()(const Parameter& t) const {
		return this->ref_ + t * this->direction_;
	}

private:
	vector_type ref_;
	direction_type direction_;

};

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
template<typename T, std::size_t Dimension>
class segment: geometry<segment_tag> {
public:
	typedef T value_type;
	typedef typename vector_type<T, Dimension>::type vector_type;
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

	template<typename Vector>
	segment(const Vector& start, const Vector& end) :
			edge_() {
		assign(start, this->edge_[GK::StartEdge]);
		assign(end, this->edge_[GK::EndEdge]);
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

	const vector_type& operator[](std::size_t index) const {
		return this->edge_[index];
	}

	vector_type& operator[](std::size_t index) {
		return this->edge_[index];
	}

	template<typename Parameter>
	vector_type operator()(const Parameter& t) const {
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

template<typename T, std::size_t Dimension>
direction<Dimension> direction_of(const line<T, Dimension>& l) {
	return l.direction();
}

template<typename T, std::size_t Dimension>
direction<Dimension> direction_of(const segment<T, Dimension>& l) {
	const typename vector_type<T, Dimension>::type v = l[GK::EndEdge]
			- l[GK::StartEdge];
	return direction_of(v);
}

template<typename T, std::size_t Dimension>
T length(const segment<T, Dimension>& l) {
	const typename vector_type<T, Dimension>::type v = l.start() - l.end();
	return norm(v);
}

template<typename T, std::size_t Dimension>
std::pair<T, T> domain(const segment<T, Dimension>&) {
	return std::make_pair(T(GK_FLOAT_ZERO), T(GK_FLOAT_ONE));
}

template<typename Vector>
aabb<Vector> boundary(const segment<Vector>& x) {
	return aabb<Vector>(x.start(), x.end());
}

template<typename Vector>
class xpolyline/*: public curve<curve_tag, Vector>*/{
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
	xpolyline() :
			X_() {

	}

	xpolyline(const xpolyline& other) :
			X_(other.X_) {

	}

	template<typename InputIterator>
	polyline(InputIterator first, InputIterator last) :
	X_(first, last) {

	}

	~xpolyline() {

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
	xpolyline subdivide(const Parameter& t, gkselection select = GK::Upper) {

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

	xpolyline& operator=(const xpolyline& rhs) {
		if (&rhs == this) {
			return *this;
		}

		this->X_ = rhs.X_;
		return *this;
	}

private:
	container_type X_;
};

}
// namespace gk

#endif /* PRIMITIVE_LINE_H_ */
