// OptimiserNonParHen.h
//
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "OptimiserVariableMetric.h"
#include "optimisable.h"

#include <limits>
#include <armadillo>

//////////////////////////////////////////////////

class OPTIMISER COptimiserNonParHen : public BOptimiserVariableMetric
{
public:
	COptimiserNonParHen () {}
	~COptimiserNonParHen () {}

	void SetParameters (const std::vector<variant>& /*parameters*/) override {};

	arma::mat UpdateHessian () override { return arma::mat (); };

	void Optimise (double target) override;

private:

};

//////////////////////////////////////////////////
