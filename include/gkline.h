/*
 * gkline.h
 *
 *  Created on: 2016/01/24
 *      Author: tmakimoto
 */

#ifndef INCLUDE_GKLINE_H_
#define INCLUDE_GKLINE_H_

namespace gk {

//template<typename Line>
//typename vector_traits<typename geometry_traits<Line>::vector_type>::direction direction_of(
//		const Line&);

template<typename Line>
void set_line(Line& l,
		const typename geometry_traits<Line>::vector_type& reference,
		const direction<geometry_traits<Line>::Dimension>& u) {
}

}  // namespace gk

#endif /* INCLUDE_GKLINE_H_ */
