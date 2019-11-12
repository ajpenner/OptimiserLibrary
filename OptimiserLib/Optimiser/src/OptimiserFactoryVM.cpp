// OptimiserFactory.h
//
// 2017
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif
#include "OptimiserFactoryVM.h"

//////////////////////////////////////////////////

OPTIMISER std::unique_ptr<IOptimiserVariableMetric> OptimiserFactoryVM::make_optimiser (eOptimiserVM::type choice)
{
	switch (choice)
	{
	case eOptimiserVM::eHennig:
		return std::make_unique<COptimiserHennig> ();
	case eOptimiserVM::eNonParHen:
		return std::make_unique<COptimiserNonParHen> ();
	case eOptimiserVM::eVMBFGS:
		return std::make_unique<COptimiserVMBFGS> ();
	case eOptimiserVM::eVMBroyden:
		return std::make_unique<COptimiserVMBroyden> ();
	case eOptimiserVM::eVMBlend:
		return std::make_unique<COptimiserVMBlend> ();
	case eOptimiserVM::eVMDFP:
		return std::make_unique<COptimiserVMDFP> ();
	case eOptimiserVM::eVMSR1:
		return std::make_unique<COptimiserVMSR1> ();
//	case eOptimiserVM::eLBFGS:
//		return std::make_unique<COptimiserLBFGS> (); // no source code??
	case eOptimiserVM::eVMPSB:
		return std::make_unique<COptimiserVMPSB> ();
	default:
		return nullptr;
	}
}

//////////////////////////////////////////////////
