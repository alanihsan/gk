/*
 * hyperplane.h
 *
 *  Created on: 2016/01/26
 *      Author: makitaku
 */

#ifndef PRIMITIVE_HYPERPLANE_H_
#define PRIMITIVE_HYPERPLANE_H_

namespace gk {

template<typename Vector>
class hyperplane {
public:
	typedef Vector vector_type;

	static const size_t Dimension = vector_traits<Vector>::Dimension;

public:
	hyperplane() {
	}

private:
	vector_type reference_;
	direction<Dimension> normal_;
};

}  // namespace gk

#endif /* PRIMITIVE_HYPERPLANE_H_ */
