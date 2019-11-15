// STLExtension.h
//
// 2016
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif

#include <cassert>
#include <algorithm>
#include <boost/range.hpp>
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

template<typename... Containers>
auto zip(Containers&&... containers) ->
	boost::iterator_range<boost::zip_iterator<decltype(boost::make_tuple(boost::make_tuple(std::begin(containers)...)))>>
{
	auto itBegin = boost::make_zip_iterator( boost::make_tuple( std::begin(containers)... ));
	auto itEnd   = boost::make_zip_iterator( boost::make_tuple( std::end(containers)... ));
	return boost::make_iterator_range(itBegin, itEnd);
}

//////////////////////////////////////////////////

// Not a proud momement
// Handles a maximum of 9 objects

#ifndef NUM_ARGS
#define NUM_ARGS0(X,_20,_19,_18,_17,_16,_15,_14,_13,_12,_11,_10,_9,_8,_7,_6,_5,_4,_3,_2,_1, N, ...) N
#define NUM_ARGS(...) NUM_ARGS0(0, __VA_ARGS__, 20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0)
#endif

#define UNZIP_1(N,tuple,x)     auto x = boost::get <N-1> (tuple)
#define UNZIP_2(N,tuple,x,...) auto x = boost::get <N-2> (tuple); UNZIP_1(N,tuple,__VA_ARGS__)
#define UNZIP_3(N,tuple,x,...) auto x = boost::get <N-3> (tuple); UNZIP_2(N,tuple,__VA_ARGS__)
#define UNZIP_4(N,tuple,x,...) auto x = boost::get <N-4> (tuple); UNZIP_3(N,tuple,__VA_ARGS__)
#define UNZIP_5(N,tuple,x,...) auto x = boost::get <N-5> (tuple); UNZIP_4(N,tuple,__VA_ARGS__)
#define UNZIP_6(N,tuple,x,...) auto x = boost::get <N-6> (tuple); UNZIP_5(N,tuple,__VA_ARGS__)
#define UNZIP_7(N,tuple,x,...) auto x = boost::get <N-7> (tuple); UNZIP_6(N,tuple,__VA_ARGS__)
#define UNZIP_8(N,tuple,x,...) auto x = boost::get <N-8> (tuple); UNZIP_7(N,tuple,__VA_ARGS__)
#define UNZIP_9(N,tuple,x,...) auto x = boost::get <N-9> (tuple); UNZIP_8(N,tuple,__VA_ARGS__)

#define UNZIP1(N,tuple,...) UNZIP_ ## N(N,tuple,__VA_ARGS__)
#define UNZIP0(N,tuple,...) UNZIP1(N,tuple,__VA_ARGS__)
#define unzip(tuple,...) UNZIP0(NUM_ARGS(__VA_ARGS__),tuple,__VA_ARGS__)

//////////////////////////////////////////////////

// The first argument is the one to be replaced by the ZipFunc results
template<typename Modified, typename T, typename U, typename ZipFunc>
void ZippedTransform (Modified& container1, T& container2, U& container3, ZipFunc zip_func)
{
	for (auto zi : zip(container1, container2, container3))
	{
		unzip( zi, &modifiedContainer, container1, container2 );
		modifiedContainer = zip_func(container1, container2);
	}
}

//////////////////////////////////////////////////

// The first argument is the one to be replaced by the ZipFunc results
template< typename Modified, typename T, typename U, typename V, typename ZipFunc>
void ZippedTransform (Modified& container1, T& container2, U& container3, V& container4, ZipFunc zip_func)
{
	for (auto zi : zip(container1, container2, container3, container4))
	{
		unzip( zi, &modifiedContainer, container1, container2, container3 );
		modifiedContainer = zip_func(container1, container2, container3);
	}
}

//////////////////////////////////////////////////

template <typename U, typename T, typename Compare>
std::vector<U> sort_permutation (
	const T& vec,
	const Compare& compare)
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
	std::transform (p.begin(), p.end(), sorted_vec.begin(), [&](U i) { return vec[i]; });
	return sorted_vec;
}

//////////////////////////////////////////////////
