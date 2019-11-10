// LineSearchFactory.h
//
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "LineSearch.h"
#include "StrongWolfe.h"
#include "WeakWolfe.h"
#include "BackTracking.h"

#include <memory>

//////////////////////////////////////////////////

namespace eLineSearch
{
	enum type
	{
		eStrongWolfe = 0,
		eWeakWolfe,
		eBackTracking,
		eMax = 255
	};
}

//////////////////////////////////////////////////

class LineSearchFactory
{
public:
	OPTIMISER static std::unique_ptr<ILineSearch> make_lineSearch (eLineSearch::type choice);
};

//////////////////////////////////////////////////