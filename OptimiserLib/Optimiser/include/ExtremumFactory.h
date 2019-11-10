// ExtremumFactory.h
//
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "MaximiserBrent.h"
#include "MinimiserBisection.h"
#include "MinimiserBoxed.h"
#include "MinimiserBrent.h"
#include "MinimiserExhaustive.h"
#include "MinimiserGolden.h"

#include <memory>

//////////////////////////////////////////////////

namespace eExtremum
{
	enum type
	{
		eMaxBrent = 0,
		eMinBrent,
		eBisection,
		eBoxed,
		eExhaustive,
		eGolden,
		eMax = 255
	};
}

//////////////////////////////////////////////////

class ExtremumFactory
{
public:
	OPTIMISER static std::unique_ptr<BExtremum> make_extremum (eExtremum::type choice);
};

//////////////////////////////////////////////////