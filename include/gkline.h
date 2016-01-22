/*
 * gkline.h
 *
 *  Created on: 2016/01/22
 *      Author: makitaku
 */

#ifndef GKLINE_H_
#define GKLINE_H_

#include "gkdirection.h"

namespace gk {

template<typename Line>
direction<geometry_traits<Line>::Dimension> direction_of(const Line& l);

}  // namespace gk

#endif /* GKLINE_H_ */
