// OptimiserVMPSB.h
//
// 2016
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "OptimiserVariableMetric.h"
#include "optimisable.h"

#include <limits>
#include <armadillo>

//////////////////////////////////////////////////

class CLinearFunction;
class CMinimiserGolden;

//////////////////////////////////////////////////

class OPTIMISER COptimiserVMPSB : public BOptimiserVariableMetric
{
public:
	COptimiserVMPSB() {}
	~COptimiserVMPSB() {}

	void SetParameters (const std::vector<variant>& parameters) override;

	arma::mat UpdateHessian () override;

	void Optimise(double target) override;

	double GetMinimalScaleValue (std::shared_ptr<CLinearFunction> vf, CMinimiserGolden& minimiser, const arma::vec& gradient);

private:

};

//////////////////////////////////////////////////
