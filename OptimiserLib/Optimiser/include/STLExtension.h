// OptimiserGRG.cpp
//
// 2016
//////////////////////////////////////////////////

#include "stdafx.h"
#include <algorithm>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/numeric/conversion/cast.hpp>

//////////////////////////////////////////////////

template<typename T>
void pop_front (std::vector<T>& vec)
{
	assert (!vec.empty ());
	vec.erase (vec.begin ());
}

//////////////////////////////////////////////////

template<typename T, typename U, typename V, typename ZipFunc>
void ZippedTransform (T& container1, U& container2, V& container3, ZipFunc zip_func)
{
	auto itBegin = boost::make_zip_iterator (boost::make_tuple (container1.begin (), container2.begin (), container3.begin ()));
	auto itEnd = boost::make_zip_iterator (boost::make_tuple (container1.end (), container2.end (), container3.end ()));
	enum { eContainer1, eContainer2, eContainer3, };

	for (auto it = itBegin; it != itEnd; ++it)
	{
		it->get<eContainer3> () = zip_func (it->get<eContainer1> (), it->get<eContainer2> ()> ());
	}

	//	std::transform(	itBegin, itEnd, itBegin, zip_func() );
}

//////////////////////////////////////////////////

template< typename T, typename U, typename V, typename W, typename ZipFunc>
void ZippedTransform (T& container1, U& container2, V& container3, W& container4, ZipFunc zip_func)
{
	auto itBegin = boost::make_zip_iterator (boost::make_tuple (container1.begin (), container2.begin (), container3.begin (), container4.begin ()));
	auto itEnd = boost::make_zip_iterator (boost::make_tuple (container1.end (), container2.end (), container3.end (), container4.end ()));
	enum { eContainer1, eContainer2, eContainer3, eContainer4, };

	for (auto it = itBegin; it != itEnd; ++it)
	{
		it->get<eContainer4>() = zip_func(it->get<eContainer1>(), it->get<eContainer2>(), it->get<eContainer3>());
	}

//	std::transform( itBegin, itEnd, itBegin, zip_func);
}

//////////////////////////////////////////////////

template <typename U, typename T, typename Compare>
std::vector<U> sort_permutation (
	const T& vec,
	Compare& compare)
{
	std::vector<U> p(vec.size());
	std::iota(p.begin (), p.end (), 0);
	std::sort (p.begin (), p.end (),
		[&](U i, U j) { return compare (vec[i], vec[j]); });
	return p;
}

//////////////////////////////////////////////////

template <typename U, typename T>
T apply_permutation (
	const T& vec,
	const std::vector<U>& p)
{
	T sorted_vec(boost::numeric_cast<U>(p.size()));
	assert (vec.size() == p.size());
	std::transform (p.begin(), p.end(), stdext::unchecked_array_iterator<decltype(sorted_vec.begin())>(sorted_vec.begin()), [&](U i) { return vec[i]; });
	return sorted_vec;
}

//////////////////////////////////////////////////
