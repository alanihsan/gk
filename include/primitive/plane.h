/*
 * plane.h
 *
 *  Created on: 2014/03/25
 *      Author: Takuya Makimoto
 */

#ifndef PRIMITIVE_PLANE_H_
#define PRIMITIVE_PLANE_H_

#include "../gkvector.h"
#include "../gkgeometry.h"

namespace gk {

/**
 * @brief This class represents plane.@n
 * @f{eqnarray*}{
 *  \mathbf{n} \cdot (\mathbf{r}(s, t) - \mathbf{r}_{ref}) = 0 \\
 *  \mathbf{r}(s, t) = \mathbf{r}_{ref} + s \mathbf{u} + t \mathbf{v}
 *  @f}
 *
 * @author Takuya Makimoto
 * @date 2015
 */
template<typename Vector>
class plane: public geometry<plane_tag, Vector> {
public:
	typedef Vector vector_type;
	typedef typename vector_traits<Vector>::value_type value_type;

	static const size_t Dimension = vector_traits<Vector>::Dimension;

	typedef direction<Dimension> direction_type;

public:

	/**
	 * @brief The default constructor.
	 */
	plane() :
			ref_(), normal_() {
	}

	/**
	 * @brief The copy constructor.
	 * @param other
	 */
	plane(const plane& other) :
			ref_(other.ref_), normal_(other.normal_) {
	}

	/**
	 *
	 * @param reference
	 * @param normal
	 */
	plane(const vector_type& reference, const direction_type& normal) :
			ref_(reference), normal_(normal) {

	}

//	plane(const vector_type& reference, const direction_type& u,
//			const direction_type& v) :
//			ref_(reference), normal_(cross(u, v)) {
//	}

	/**
	 * @brief Creates from 3 position vectors.
	 *
	 * @param first
	 * @param second
	 * @param third
	 */
	plane(const vector_type& first, const vector_type& second,
			const vector_type& third) :
			ref_(first), normal_(normal(second - first, third - first)) {
	}

	/**
	 * @brief Destructor.
	 */
	~plane() {
	}

	/**
	 * @brief Returns a constant reference to the reference point.
	 * @return
	 */
	const vector_type& reference() const {
		return this->ref_;
	}

	/**
	 * @brief Returns a mutable reference to the reference point.
	 * @return
	 */
	vector_type& reference() {
		return this->ref_;
	}

	/**
	 * @brief Returns a constant reference to the normal vector.
	 * @return
	 */
	const direction_type& normal() const {
		return this->normal_;
	}

	direction_type& normal() {
		return this->normal_;
	}

	/**
	 * @brief Computes the @f$u@f$ direction.
	 * @return
	 */
	direction_type u_direction() const {
		return this->u_direction_();
	}

	/**
	 * @brief Computes the @f$v@f$ direction.
	 * @return
	 */
	direction_type v_direction() const {
		return this->v_direction_();
	}

	/**
	 * @brief Computes a position at @f$(s, t)@f$.
	 *
	 * @tparam Parameter Type of parameter.
	 * @param s u-coordinate on this plane. (Local coordinates)
	 * @param t v-coordinate on this plane. (Local coordinates)
	 * @return Position vector in the global system.
	 */
	template<typename Parameter>
	vector_type operator()(const Parameter& s, const Parameter& t) const {
		return this->ref_ + s * this->u_direction_() + t * this->v_direction_();
	}

	plane& operator=(const plane& rhs) {
		if (&rhs == this) {
			return *this;
		}

		this->ref_ = rhs.ref_;
		this->normal_ = rhs.normal_;
		return *this;
	}

private:
	vector_type ref_; ///< The reference position of this plane.
	direction_type normal_; ///< The normal vector of this plane.

private:
	direction_type u_direction_() const {
		if (this->normal_[GK::X]
				== typename direction_type::value_type(GK_FLOAT_ZERO)) {
			return basis<Dimension>()[GK::X];
		}

		if (this->normal_[GK::Y]
				== typename direction_type::value_type(GK_FLOAT_ZERO)) {
			return basis<Dimension>()[GK::Y];
		}

		if (this->normal_[GK::Z]
				== typename direction_type::value_type(GK_FLOAT_ZERO)) {
			return basis<Dimension>()[GK::Z];
		}

		vector_type u;
		u[GK::X] = value_type(-this->normal_[GK::Y]);
		u[GK::Y] = value_type(this->normal_[GK::X]);
		u[GK::Z] = value_type(GK_FLOAT_ZERO);
		return direction_type(u);
	}

	direction_type v_direction_() const {
		const value_type unit = value_type(GK_FLOAT_ONE);
		return normal_direction(unit * this->normal_,
				unit * this->u_direction_());
	}
};

template<typename Vector>
direction<vector_traits<Vector>::Dimension> direction_of(
		const plane<Vector>& x) {
	return x.normal();
}

} // namespace gk

#endif /* PRIMITIVE_PLANE_H_ */
