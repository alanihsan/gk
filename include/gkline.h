/*
 * gkline.h
 *
 *  Created on: 2016/01/24
 *      Author: tmakimoto
 */

#ifndef INCLUDE_GKLINE_H_
#define INCLUDE_GKLINE_H_

namespace gk {

template<typename Line>
typename vector_traits<typename geometry_traits<Line>::vector_type>::direction direction_of(
		const Line&);

}  // namespace gk

#endif /* INCLUDE_GKLINE_H_ */
