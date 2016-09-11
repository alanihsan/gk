/*
 * gkdef.h
 *
 *  Created on: 2013/04/09
 *      Author: Takuya Makimoto
 */

/**
 * \file
 * \brief
 */

#ifndef GKDEF_H_
#define GKDEF_H_

#include <cstdlib>
#include <stdint.h>
#include <limits>

#include "config/gkconfig.h"

///*
// * Integer Value Type．
// */
//#ifdef GK_SIZEOF_INT
//#	if GK_SIZEOF_INT == __SIZEOF_INT__
//#		define GK_INT_TYPE  int32_t
//#		define GK_UINT_TYPE uint32_t
//#		define GK_INT_MAX   __INT32_MAX__
//#		define GK_UINT_MAX  __UINT32_MAX__
//# 	elif GK_SIZEOF_INT == __SIZEOF_LONG_LONG__
//#		define GK_INT_TYPE  int64_t
//#		define GK_UINT_TYPE uint64_t
//#		define GK_INT_MAX   __INT64_MAX__
//#		define GK_UINT_MAX  __UINT64_MAX__
//#	else
//#		error "GK_INT_BYTE" must set 4 or 8.
//#	endif
//#else
//#	define GK_SIZEOF_INT 4
//#	define GK_INT_TYPE int32_t
//#	define GK_UINT_TYPE uint32_t
//#	define GK_INT_MAX   __INT32_MAX__
//#	define GK_UINT_MAX  __UINT32_MAX__
//#endif
//typedef GK_INT_TYPE gkint; ///< Signed integer.
//typedef GK_UINT_TYPE gkuint; ///< Unsigned integer.
//#undef GK_INT_TYPE
//#undef GK_UINT_TYPE
//
//typedef gkuint gksize;

typedef bool gkselection;

typedef bool gkminmax;
const bool GK_Min = false;
const bool GK_Max = true;

#define GK_FLOAT_ZERO 0.0
#define GK_FLOAT_ONE 1.0
#define GK_FLOAT_NEGATIVE_ONE -1.0

namespace gk {
/*
 * Floating Value Type．
 */
#ifdef GK_SIZEOF_FLOAT
#	if GK_SIZEOF_FLOAT == __SIZEOF_FLOAT__
#		define GK_FLOAT_TYPE float
#	elif GK_SIZEOF_FLOAT == __SIZEOF_DOUBLE__
#		define GK_FLOAT_TYPE double
#	else
#		error "GK_FLOAT_TYPE" must set 4 or 8.
#	endif
#else
#	define GK_SIZEOF_FLOAT 8
#		define GK_FLOAT_TYPE double
#endif
typedef GK_FLOAT_TYPE float_type;
#undef GK_FLOAT_TYPE

template<std::size_t N>
struct number_tag {
	static const std::size_t Value = N;
};

template<std::size_t Dimension>
struct dimension_tag: public number_tag<Dimension> {
	static const std::size_t Size = Dimension;
};

template<typename Scalar>
struct scalar_traits {
	static const std::size_t Precision = 0;
};

template<std::size_t Precision>
struct float_traits {
	typedef void value_type;
};

template<>
struct float_traits<4> {
	typedef float value_type;
};

template<>
struct float_traits<8> {
	typedef double value_type;
};

template<>
struct float_traits<16> {
	typedef long double value_type;
};

/**
 * @brief Struct of constant objects.
 */
struct GK {

//#ifdef GK_EIGEN_ENABLED
//	static const gksize DynamicSize = static_cast<gksize>(Eigen::Dynamic);
//#endif

	/**
	 * @brief Enum about space dimensions.
	 */
	typedef enum {
		GK_NonDim = 0, ///< Non-dimension.
		GK_1D = 1, ///< 1 dimension.
		GK_2D = 2, ///< 2 dimension.
		GK_3D = 3 ///< 3 dimension.
	} DimensionNumber;

	static const std::size_t D3 = 3;

	typedef enum {
		GK_SinglePrecision = __SIZEOF_FLOAT__,
		GK_DoublePrecision = __SIZEOF_DOUBLE__,
		GK_LongDoublePrecision = __SIZEOF_LONG_DOUBLE__,
	} Precision;

	typedef enum {
		StartEdge, EndEdge, EdgeSize
	} Edge;

	static const std::size_t X = 0;
	static const std::size_t Y = 1;
	static const std::size_t Z = 2;
	static const std::size_t W = 3;

	static const std::size_t U = 0;
	static const std::size_t V = 1;

	static const std::size_t MinTag = 0;
	static const std::size_t MaxTag = 1;

//	static const gkfloat Nan;
//	static const gkfloat Zero;
//	static const gkfloat Unit;
//	static const gkfloat Infinity;

	static const bool Upper = false;
	static const bool Lower = true;

	static const bool InfiniteLength = false;
	static const bool FiniteLength = true;
};

}  // namespace gk

#endif /* GKDEF_H_ */
