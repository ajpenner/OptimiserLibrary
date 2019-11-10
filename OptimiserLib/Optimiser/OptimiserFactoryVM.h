// OptimiserFactory.h
//
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "OptimiserVariableMetric.h"
#include "OptimiserHennig.h"
#include "OptimiserVMLBFGS.h"
#include "OptimiserNonParHen.h"
#include "OptimiserVMBFGS.h"
#include "OptimiserVMBroyden.h"
#include "OptimiserVMBlend.h"
#include "OptimiserVMDFP.h"
#include "OptimiserVMSR1.h"
#include "OptimiserVMPSB.h"

#include <memory>

//////////////////////////////////////////////////

namespace eOptimiserVM
{
	enum type
	{
		eHennig,
		eNonParHen,
		eVMBFGS,
		eVMBroyden,
		eVMBlend,
		eVMDFP,
		eVMSR1,
		eVMLBFGS, // not implemented
		eVMPSB,
		eMax = 255
	};
}

//////////////////////////////////////////////////

class OptimiserFactoryVM
{
public:
	OPTIMISER static std::unique_ptr<IOptimiserVariableMetric> make_optimiser (eOptimiserVM::type choice);
};

//////////////////////////////////////////////////