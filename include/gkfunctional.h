/*
 * gkfunctional.h
 *
 *  Created on: 2016/01/18
 *      Author: makitaku
 */

#ifndef GKFUNCTIONAL_H_
#define GKFUNCTIONAL_H_

namespace gk {

template<typename T1, typename T2>
struct multiplies_result {
	typedef void value_type;
};

template<>
struct multiplies_result<float, float> {
	typedef float value_type;
};

template<>
struct multiplies_result<double, double> {
	typedef double value_type;
};

template<>
struct multiplies_result<long double, long double> {
	typedef long double value_type;
};

template<>
struct multiplies_result<float, double> {
	typedef float value_type;
};

template<>
struct multiplies_result<double, float> {
	typedef float value_type;
};

template<>
struct multiplies_result<float, long double> {
	typedef float value_type;
};

template<>
struct multiplies_result<double, long double> {
	typedef double value_type;
};

template<>
struct multiplies_result<long double, float> {
	typedef float value_typ;
};

template<>
struct multiplies_result<long double, double> {
	typedef double value_type;
};

template<typename T1, typename T2>
struct divides_result {
	typedef void value_type;
};

template<>
struct divides_result<float, float> {
	typedef float value_type;
};

template<>
struct divides_result<double, double> {
	typedef double value_type;
};

template<>
struct divides_result<long double, long double> {
	typedef long double value_type;
};

template<>
struct divides_result<float, double> {
	typedef float value_type;
};

template<>
struct divides_result<double, float> {
	typedef float value_type;
};

template<>
struct divides_result<float, long double> {
	typedef float value_type;
};

template<>
struct divides_result<double, long double> {
	typedef double value_type;
};

template<>
struct divides_result<long double, float> {
	typedef float value_type;
};

template<>
struct divides_result<long double, double> {
	typedef double value_type;
};

template<typename T1, typename T2, typename Result = typename multiplies_result<
		T1, T2>::value_type> struct multiplies {
	typedef T1 first_argument_type;
	typedef T2 second_argument_type;
	typedef Result result_type;

	Result operator()(const T1& x, const T2& y) const {
		return x * y;
	}
};

template<typename T1, typename T2, typename Result = typename divides_result<T1,
		T2>::value_type> struct divides {
	typedef T1 first_argument_type;
	typedef T2 second_argument_type;
	typedef Result result_type;

	Result operator()(const T1& x, const T2& y) const {
		return x / y;
	}
};

}  // namespace gk

#endif /* GKFUNCTIONAL_H_ */
