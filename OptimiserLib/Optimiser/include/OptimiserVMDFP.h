// OptimiserVariableMetric.h
//
// 2016
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "optimiserVariableMetric.h"
#include "optimisable.h"

#include <limits>
#include <armadillo>

//////////////////////////////////////////////////

class CLinearFunction;
class CMinimiserGolden;

//////////////////////////////////////////////////

class OPTIMISER COptimiserVMDFP : public BOptimiserVariableMetric
{
public:
	COptimiserVMDFP () {}
	~COptimiserVMDFP () {}

	void SetParameters (const std::vector<variant>& /*parameters*/) override {};

	void Optimise(double target) override;

	arma::mat UpdateHessian() override;

	double GetMinimalScaleValue (std::shared_ptr<CLinearFunction> vf, CMinimiserGolden& minimiser, const arma::vec& gradient);

private:

};

//////////////////////////////////////////////////