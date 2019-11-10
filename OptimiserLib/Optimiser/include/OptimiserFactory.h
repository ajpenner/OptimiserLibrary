// OptimiserFactory.h
//
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

// Constrained
#include "OptimiserEllipsoid.h"
#include "OptimiserGRG.h"
#include "OptimiserPenalty.h"
#include "OptimiserMOM.h"
#include "OptimiserComplex.h"
#include "OptimiserRandom.h"

// Unconstrained
#include "OptimiserNewton.h"
#include "OptimiserBrent.h"
#include "OptimiserPowell.h"
#include "OptimiserBox.h"
#include "OptimiserSimplex.h"
#include "OptimiserHJ.h"
#include "OptimiserCauchy.h"
#include "OptimiserMarquardt.h"
#include "OptimiserCG.h"

#include <memory>

//////////////////////////////////////////////////

namespace eOptimiser
{
	enum type
	{
		// Derivative
		eNewton,
		eCauchy,
		eMarquardt,
		eConjugateGradient,
		// No derivative
		eBrent,
		ePowell,
		eBox,
		eSimplex,
		eHJ,
		eMax = 255
	};
}

//////////////////////////////////////////////////

namespace eConstraintOptimiser
{
	enum type
	{
		eEllipsoid = 0,
		eGRG,
		ePenalty,
		ePenaltyMOM,
		eComplex, // does not work
		eRandom,
		eMax = 255
	};
}

//////////////////////////////////////////////////

class OptimiserFactory
{
public:
	OPTIMISER static std::unique_ptr<IOptimiser> make_optimiser (eOptimiser::type choice);
};

//////////////////////////////////////////////////

class ConstraintOptimiserFactory
{
public:
	OPTIMISER static std::unique_ptr<BConstraintOptimiser> make_optimiser (eConstraintOptimiser::type choice);
};

//////////////////////////////////////////////////