/*
 * hyperplane.h
 *
 *  Created on: 2016/01/26
 *      Author: Takuya Makimoto
 */

#ifndef PRIMITIVE_HYPERPLANE_H_
#define PRIMITIVE_HYPERPLANE_H_

#include "../gkvector.h"

namespace gk {

/**
 * @brief
 *
 * @date 2016/06/28
 * @author Takuya Makimoto
 */
template<typename Vector>
class hyperplane {
public:
	typedef Vector vector_type;

	static const size_t Dimension = vector_traits<Vector>::Dimension;

	typedef direction<Dimension> direction_type;

public:
	hyperplane() {
	}

	hyperplane(const hyperplane& other) :
			reference_(other.reference_), normal_(other.normal_) {
	}

	hyperplane(const vector_type& reference, const direction_type& normal) :
			reference_(reference), normal_(normal) {
	}

	~hyperplane() {
	}

	vector_type reference() const {
		return this->reference_;
	}

	void reference(const vector_type r) {
		this->reference_ = r;
	}

	direction_type normal() const {
		return this->normal_;
	}

	void normal(const direction_type& n) {
		this->normal_ = n;
	}

	hyperplane& operator=(const hyperplane& rhs) {
		if (&rhs == this) {
			return *this;
		}

		this->reference_ = rhs.reference_;
		this->normal_ = rhs.normal_;

		return *this;
	}

	/**
	 * @brief
	 *
	 * @tparam T Type of line in 2D, type of plane in 3D.
	 *
	 * @param rhs
	 * @return
	 */
	template<typename T>
	hyperplane& operator=(const T& rhs) {
		return this->assign_(rhs, typename geometry_traits<T>::category(),
				dimension_tag<geometry_traits<T>::Dimension>());
	}

private:
	vector_type reference_;
	direction<Dimension> normal_;

private:
	template<typename Line>
	hyperplane& assign_(const Line& line, line_tag, dimension_tag<GK::GK_2D>) {
		this->reference_ = line(GK_FLOAT_ZERO);

		const direction_type t = direction_of(line);
		vector_type v;
		v[GK::X] = -t[GK::Y];
		v[GK::Y] = t[GK::X];
		const direction_type n = normalize(v);
		this->normal_ = n;

		return *this;
	}
};

}  // namespace gk

#endif /* PRIMITIVE_HYPERPLANE_H_ */
