/*
 * gkxaabb.h
 *
 *  Created on: 2016/09/01
 *      Author: tmakimoto
 */

#ifndef INCLUDE_GKXAABB_H_
#define INCLUDE_GKXAABB_H_

//#include "vector/algebra.h"
#include "xgkvector.h"

namespace gk {

/**
 * @brief Axis-aligned bounding box (AABB).
 *
 * @author Takuya Makimoto
 * @date 2014
 */
template<typename Vector, std::size_t DimensionSize, typename T>
class xaabb {
public:
	typedef typename vector_type<DimensionSize, T>::type vector_type;
//	typedef Eigen::Matrix<T,1,DimensionSize,Eigen::RowMajor> vector_type;

	static const std::size_t Dimension = DimensionSize;
public:
	/**
	 * @brief Default constructor.
	 */
	xaabb() {
	}

	/**
	 * @brief Copy constructor.
	 * @param other
	 */
	xaabb(const xaabb& other) :
			min_(other.min_), max_(other.max_) {
	}

	xaabb(const vector_type& u, const vector_type& v) :
			min_(), max_() {
		this->set_(u, v, dimension_tag<Dimension>());
	}

	template<class InputIterator>
	xaabb(InputIterator first, InputIterator last) :
			min_(), max_() {
		typedef InputIterator iterator;
//		typedef std::iterator_traits<InputIterator> traits;

		{
			iterator p = first;
			this->min_ = *p;
			this->max_ = *p;

			++p;

			while (p != last) {
				for (size_t i = 0; i < Dimension; ++i) {
					this->min_[i] = std::min(this->min_[i], (*p)[i]);
					this->max_[i] = std::max(this->max_[i], (*p)[i]);
				}
				++p;
			}
		}
	}

	/**
	 * @brief Destructor.
	 */
	~xaabb() {
	}

	const vector_type& min() const {
		return this->min_;
	}

	const vector_type& max() const {
		return this->max_;
	}

	xaabb& operator=(const xaabb& rhs) {
		if (&rhs == this) {
			return *this;
		}

		this->min_ = rhs.min_;
		this->max_ = rhs.max_;

		return *this;
	}

private:
	vector_type min_;
	vector_type max_;

private:

	template<std::size_t Dimension>
	void set_(const vector_type& u, const vector_type& v,
			dimension_tag<Dimension>) {
		for (std::size_t i = 0; i < Dimension; ++i) {
			(u[i] < v[i]) ?
					this->min_[i] = u[i], this->max_[i] = v[i] : this->min_[i] =
							v[i], this->max_[i] = u[i];
		}
	}

	void set_(const vector_type& u, const vector_type& v,
			dimension_tag<GK::GK_2D>) {
		(u[GK::X] < v[GK::X]) ?
				(this->min_[GK::X] = u[GK::X], this->max_[GK::X] = v[GK::X]) :
				(this->min_[GK::X] = v[GK::X], this->max_[GK::X] = u[GK::X]);

		(u[GK::Y] < v[GK::Y]) ?
				(this->min_[GK::Y] = u[GK::Y], this->max_[GK::Y] = v[GK::Y]) :
				(this->min_[GK::Y] = v[GK::Y], this->max_[GK::Y] = u[GK::Y]);
	}

	void set_(const vector_type& u, const vector_type& v,
			dimension_tag<GK::GK_3D>) {
		set_(u, v, dimension_tag<GK::GK_2D>());

		(u[GK::Z] < v[GK::Z]) ?
				this->min_[GK::Z] = u[GK::Z], this->max_[GK::Z] = v[GK::Z] :
				this->min_[GK::Z] = v[GK::Z], this->max_[GK::Z] = u[GK::Z];
	}
};

}  // namespace gk

#endif /* INCLUDE_GKXAABB_H_ */
