// OptimiserVMBFGS.h
//
// 2016
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "optimiser.h"
#include "optimisable.h"
#include "OptimiserVariableMetric.h"
#include <armadillo>

//////////////////////////////////////////////////

class CLinearFunction;
class CMinimiserGolden;

//////////////////////////////////////////////////

class OPTIMISER COptimiserVMSR1 : public BOptimiserVariableMetric
{
public:
	COptimiserVMSR1 () {}
	~COptimiserVMSR1 () {}

	void SetParameters (const std::vector<variant>& /*parameters*/) override {};

	void Optimise(double target) override;

	arma::mat UpdateHessian () override;

	double GetMinimalScaleValue (std::shared_ptr<CLinearFunction> ptrVF, CMinimiserGolden& minimiser, const arma::vec& direction);

private:

};

//////////////////////////////////////////////////