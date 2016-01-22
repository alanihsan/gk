/*
 * plane.h
 *
 *  Created on: 2014/03/25
 *      Author: Takuya Makimoto
 */

#ifndef PRIMITIVE_PLANE_H_
#define PRIMITIVE_PLANE_H_

#include "../gkvector.h"

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
class plane/*: public surface<plane_tag, Vector, Parameter>*/{
public:
	typedef Vector vector_type;
	typedef typename vector_traits<Vector>::value_type value_type;

	static const size_t Dimension = vector_traits<Vector>::Dimension;

	typedef typename vector_traits<Vector>::direction_type direction_type;

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

//	plane(const gkinvalid& invalid) :
//			ref_(invalid), normal_(invalid) {
//	}

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
			ref_(first), normal_(
					cross<direction_type, Vector, Vector>(
							normalize(second - first) / norm(second - first),
							normalize(third - first))) {
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

//	/**
//	 * @brief Computes the position projected to this plane in @a direction.
//	 * @param position
//	 * @param direction
//	 * @return
//	 */
//	gkvector<GK::GK_2D> project(const gkvector<GK::GK_3D>& position,
//			const gkdirection<GK::GK_3D>& direction) const {
//		return gkinvalid();
//	}

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

//	direction_type operator[](size_t n) const {
//
//	}

//	/**
//	 * @brief Computes the position in the global system.
//	 *
//	 * @param u u-coordinate on this plane. (Local coordinates)
//	 * @param v v-coordinate on this plane. (Local coordinates)
//	 * @return
//	 */
//	template<typename Parameter>
//	vector_type operator()(const Parameter& t) const {
//		return this->ref_ + t[GK::U] * this->u_direction_()
//				+ t[GK::V] * this->v_direction_();
//	}

	/**
	 * @brief Compute a position at @f$(s, t)@f$.
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

//	bool operator()(const vector_type& p) const {
//		return dot(p - this->ref_, this->normal_) >= parameter(GK_FLOAT_ZERO);
//	}

	plane& operator=(const plane& rhs) {
		if (&rhs == this) {
			return *this;
		}

		this->ref_ = rhs.ref_;
		this->normal_ = rhs.normal_;
		return *this;
	}

private:
	vector_type ref_;
	direction_type normal_;

private:
	direction_type u_direction_() const {
		if (this->normal_[GK::X] == direction_type::value_type(GK_FLOAT_ZERO)) {
//			return const_direction<vector_type>()[GK::X];
			return basis<Vector>()[GK::X];
		}

		if (this->normal_[GK::Y] == direction_type::value_type(GK_FLOAT_ZERO)) {
			return basis<Vector>()[GK::Y];
		}

		if (this->normal_[GK::Z] == direction_type::value_type(GK_FLOAT_ZERO)) {
			return basis<Vector>()[GK::Z];
		}

		return direction_type(
				vector_type(value_type(-this->normal_[GK::Y]),
						value_type(this->normal_[GK::X],
								value_type(GK_FLOAT_ZERO))));

//		typename direction_type::const_iterator p = std::find(
//				this->normal_.begin(), this->normal_.end(),
//				direction_type::value_type(GK_FLOAT_ZERO));
//		if (p != this->normal_.end()) {
//			return basis<Vectpr>()[std::distance(this->normal_.begin(), p)];
//		}
//
//		direction_type d = this->normal_;

	}

	direction_type v_direction_() const {
		return cross(this->normal_, this->u_direction_());
	}
};

template<typename Vector>
direction<Vector> direction_of(const plane<Vector>& x) {
	return x.normal();
}

} // namespace gk

#endif /* PRIMITIVE_PLANE_H_ */
