/*
 * gkcurve.h
 *
 *  Created on: 2015/09/25
 *      Author: makitaku
 */

#ifndef GKCURVE_H_
#define GKCURVE_H_

#include "gkgeometry.h"
#include "gkvector.h"
#include "gkaabb.h"

namespace gk {

//template<typename CurveCategory, typename Vector>
//struct curve;
//
//#define GK_CURVE_BASE_TEMPLATE_CLASS(category) \
//		template<typename Vector> \
//		struct curve<category, Vector> : public geometry<category, Vector> { \
//		};
//
//GK_CURVE_BASE_TEMPLATE_CLASS(curve_tag)
//GK_CURVE_BASE_TEMPLATE_CLASS(line_tag)
//GK_CURVE_BASE_TEMPLATE_CLASS(ray_tag)
//GK_CURVE_BASE_TEMPLATE_CLASS(segment_tag)
//GK_CURVE_BASE_TEMPLATE_CLASS(circle_tag)
//GK_CURVE_BASE_TEMPLATE_CLASS(free_curve_tag)
//#undef GK_CURVE_BASE_TEMPLATE_CLASS
//
//template<typename Curve>
//struct curve_traits: public geometry_traits<Curve> {
//	typedef typename Curve::vector_type vector_type;
//};
//
///**
// * @brief Returns true if a curve is closed, otherwise false.
// * @param curve
// * @return
// */
//template<typename Curve>
//bool is_closed(const Curve& curve);
//
///**
// // * @brief Returns a domain of a curve.
// // * @param curve
// // * @return
// // */
////template<typename Curve>
////std::pair<typename curve_traits<Curve>::parameter,
////		typename curve_traits<Curve>::parameter> domain(const Curve& curve);
///**
// * @brief
// * @param x
// * @return
// */
//template<typename Curve>
//typename curve_traits<Curve>::distance_type length(const Curve& x);
//
//template<typename Curve, typename T>
//typename divides_result<Curve, T>::value_type derivative(const Curve& x);
//
///**
// *
// * @param x A curve to be subdivided.
// * @param t A parameter.
// * @param result
// * @return
// */
//template<typename Curve, typename Parameter, typename OutputIterator>
//OutputIterator subdivide(const Curve& x, const Parameter& t,
//		OutputIterator result);

}// namespace gk

#endif /* GKCURVE_H_ */
