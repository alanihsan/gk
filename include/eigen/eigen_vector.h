/*
 * eigen_vector.h
 *
 *  Created on: 2016/01/22
 *      Author: makitaku
 */

#ifndef GKEIGEN_H_
#define GKEIGEN_H_

#include "../gkvector.h"
#include <Eigen/Core>

#define GK_EIGEN_ROWVECTOR 0
#define GK_EIGEN_COLUMNVECTOR 1

namespace gk {

template<typename Scalar, int DimensionSize, int Options>
struct vector_traits<Eigen::Matrix<Scalar, DimensionSize, 1, Options> > {
	typedef Scalar value_type;

	static const size_t Dimension = DimensionSize;
	static const bool IsHomegeneous = false;
};

template<typename Scalar, int DimensionSize, int Options>
struct vector_traits<const Eigen::Matrix<Scalar, DimensionSize, 1, Options> > {
	typedef Scalar value_type;

	static const size_t Dimension = DimensionSize;
	static const bool IsHomegeneous = false;
};

template<typename Scalar, int DimensionSize, int Options>
struct vector_traits<Eigen::Matrix<Scalar, 1, DimensionSize, Options> > {
	typedef Scalar value_type;

	static const size_t Dimension = DimensionSize;
	static const bool IsHomegeneous = false;
};

template<typename Scalar, int DimensionSize, int Options>
struct vector_traits<const Eigen::Matrix<Scalar, 1, DimensionSize, Options> > {
	typedef Scalar value_type;

	static const size_t Dimension = DimensionSize;
	static const bool IsHomegeneous = false;
};

/*
 * Traits and functions about a column vector in Eigen.
 */

template<typename Scalar1, typename Scalar2, int DimensionSize, int Options1,
		int Options2>
typename multiplies_result<Scalar1, Scalar2>::value_type dot(
		const Eigen::Matrix<Scalar1, DimensionSize, 1, Options1>& u,
		const Eigen::Matrix<Scalar2, DimensionSize, 1, Options2>& v) {
	return u.dot(v);
}

template<typename Scalar, int DimensionSize, int Options>
Scalar norm(const Eigen::Matrix<Scalar, DimensionSize, 1, Options>& v) {
	return v.norm();
}

/*
 * Traits and functions about a row vector in Eigen.
 */

template<typename Scalar1, typename Scalar2, int DimensionSize, int Options1,
		int Options2>
typename multiplies_result<Scalar1, Scalar2>::value_type dot(
		const Eigen::Matrix<Scalar1, 1, DimensionSize, Options1>& u,
		const Eigen::Matrix<Scalar2, 1, DimensionSize, Options2>& v) {
	return u.dot(v);
}

template<typename Scalar, int DimensionSize, int Options>
Scalar norm(const Eigen::Matrix<Scalar, 1, DimensionSize, Options>& v) {
	return v.norm();
}

//#if GK_USING_EIGEN==GK_EIGEN_COLUMNVECTOR
//
//template<typename Scalar>
//Eigen::Matrix<Scalar, GK::GK_2D, 1, Eigen::ColMajor> operator*(
//		const Scalar& alpha, const direction<GK::GK_2D>& u) {
//	return Eigen::Matrix<Scalar, GK::GK_2D, 1>(alpha * u[GK::X],
//			alpha * u[GK::Y]);
//}
//
//template<typename Scalar>
//Eigen::Matrix<Scalar, GK::GK_3D, 1, Eigen::ColMajor> operator*(const Scalar& alpha,
//		const direction<GK::GK_3D>& u) {
//	return Eigen::Matrix<Scalar, GK::GK_3D, 1>(alpha * u[GK::X],
//			alpha * u[GK::Y], alpha * u[GK::Z]);
//}
//
//template<typename Scalar>
//Eigen::Matrix<Scalar, GK::GK_2D, 1, Eigen::ColMajor> operator*(const direction<GK::GK_2D>& u,
//		const Scalar& alpha) {
//	return alpha * u;
//}
//
//template<typename Scalar>
//Eigen::Matrix<Scalar, GK::GK_3D, 1, Eigen::ColMajor> operator*(const direction<GK::GK_3D>& u,
//		const Scalar& alpha) {
//	return alpha * u;
//}
//
//#else

template<typename Scalar>
Eigen::Matrix<Scalar, 1, GK::GK_2D, Eigen::RowMajor> operator*(
		const Scalar& alpha, const direction<GK::GK_2D>& u) {
	return Eigen::Matrix<Scalar, 1, GK::GK_2D>(alpha * u[GK::X],
			alpha * u[GK::Y]);
}

template<typename Scalar>
Eigen::Matrix<Scalar, 1, GK::GK_3D, Eigen::RowMajor> operator*(
		const Scalar& alpha, const direction<GK::GK_3D>& u) {
	return Eigen::Matrix<Scalar, 1, GK::GK_3D>(alpha * u[GK::X],
			alpha * u[GK::Y], alpha * u[GK::Z]);
}

template<typename Scalar>
Eigen::Matrix<Scalar, 1, GK::GK_2D, Eigen::RowMajor> operator*(
		const direction<GK::GK_2D>& u, const Scalar& alpha) {
	return alpha * u;
}

template<typename Scalar>
Eigen::Matrix<Scalar, 1, GK::GK_3D, Eigen::RowMajor> operator*(
		const direction<GK::GK_3D>& u, const Scalar& alpha) {
	return alpha * u;
}
//#endif

template<typename BinaryOp, typename Vector>
struct vector_traits<Eigen::CwiseBinaryOp<BinaryOp, Vector, Vector> > {
	typedef typename vector_traits<Vector>::value_type value_type;
	static const size_t Dimension = vector_traits<Vector>::Dimension;
	static const bool IsHomegeneous = false;
};

template<typename BinaryOp, typename Vector>
struct vector_traits<const Eigen::CwiseBinaryOp<BinaryOp, Vector, Vector> > {
	typedef typename vector_traits<Vector>::value_type value_type;
	static const size_t Dimension = vector_traits<Vector>::Dimension;
	static const bool IsHomegeneous = false;
};

}  // namespace gk

#endif /* GKEIGEN_H_ */
