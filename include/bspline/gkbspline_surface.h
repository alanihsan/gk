/*
 * gkbspline_surface.h
 *
 *  Created on: 2014/09/01
 *      Author: makitaku
 */

#ifndef GKBSPLINE_SURFACE_H_
#define GKBSPLINE_SURFACE_H_

#include <bspline/knotvector.h.del>
#include <gknet.h>

class gkbspline_surface: public gkgeometry<gk_surface_tag, GK::GK_3D,
		gkvector<GK::GK_2D> > {
public:
	gkbspline_surface();
	gkbspline_surface(const gkbspline_surface& other);
	gkbspline_surface(const gkknotvector& U, const gkknotvector& V,
			const gknet& controls);
	~gkbspline_surface();

	const gkknotvector& knot_vector_u() const;
	const gkknotvector& knot_vector_v() const;
	const gknet& controls() const;

	void translate(const gkvector& translation);
	void rotate(const gkvector& rotation);
	void scale(gkfloat factor);

	gkvector<GK::GK_3D> operator()(gkparam u, gkparam v) const;
	gkvector<GK::GK_3D> operator()(const parameter& t) const;

	gkbspline_surface& operator=(const gkbspline_surface& rhs);

private:
	gkknotvector U_; /// Knot vector of @f$u@f$ direction.
	gkknotvector V_; /// Knot vector of @f$v@f$ direction.
	gknet controls_;

private:
	gkuint order_u_() const;
	gkuint order_v_() const;
};

#endif /* GKBSPLINE_SURFACE_H_ */
