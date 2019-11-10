// DirectionFactory.h
//
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "Direction.h"
#include "DirectionCauchy.h"
#include "DirectionNewton.h"
#include "DirectionQuasiNewton.h"

#include <memory>

//////////////////////////////////////////////////

namespace eDirection
{
	enum type
	{
		eCauchy = 0,
		eNewton,
		eQuasiNewton, 
		eMax = 255
	};
}

//////////////////////////////////////////////////

class DirectionFactory
{
public:
	OPTIMISER static std::unique_ptr<IDirection> make_direction (eDirection::type choice);
};

//////////////////////////////////////////////////
