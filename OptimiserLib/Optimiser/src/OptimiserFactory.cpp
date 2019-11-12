// OptimiserFactory.h
//
// 2017
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif
#include "OptimiserFactory.h"

//////////////////////////////////////////////////

OPTIMISER std::unique_ptr<IOptimiser> OptimiserFactory::make_optimiser (eOptimiser::type choice)
{
	switch (choice)
	{
	case eOptimiser::eNewton:
		return std::make_unique<COptimiserNewton> ();
	case eOptimiser::eBrent:
		return std::make_unique<COptimiserBrent> ();
	case eOptimiser::ePowell:
		return std::make_unique<COptimiserPowell> ();
	case eOptimiser::eBox:
		return std::make_unique<COptimiserBox> ();
	case eOptimiser::eSimplex:
		return std::make_unique<COptimiserSimplex> ();
	case eOptimiser::eHJ:
		return std::make_unique<COptimiserHJ> ();
	case eOptimiser::eCauchy:
		return std::make_unique<COptimiserCauchy> ();
	case eOptimiser::eMarquardt:
		return std::make_unique<COptimiserMarquardt> ();
	case eOptimiser::eConjugateGradient:
		return std::make_unique<COptimiserCG> ();
	default:
		return nullptr;
	}
}

//////////////////////////////////////////////////

OPTIMISER std::unique_ptr<BConstraintOptimiser> ConstraintOptimiserFactory::make_optimiser (eConstraintOptimiser::type choice)
{
	switch (choice)
	{
	case eConstraintOptimiser::eEllipsoid:
		return std::make_unique<COptimiserEllipsoid> ();
	case eConstraintOptimiser::eGRG:
		return std::make_unique<COptimiserGRG> ();
	case eConstraintOptimiser::ePenalty:
		return std::make_unique<COptimiserPenalty> ();
	case eConstraintOptimiser::ePenaltyMOM:
		return std::make_unique<COptimiserMOM> ();
	case eConstraintOptimiser::eComplex:
		return std::make_unique<COptimiserComplex> ();
	case eConstraintOptimiser::eRandom:
		return std::make_unique<COptimiserRandom> ();
	default:
		return nullptr;
	}
}

//////////////////////////////////////////////////
