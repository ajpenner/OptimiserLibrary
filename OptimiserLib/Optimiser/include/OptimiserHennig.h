// OptimiserHenning.h
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

class CLinearFunction;
class CMinimiserGolden;

//////////////////////////////////////////////////

class OPTIMISER COptimiserHennig : public BOptimiserVariableMetric
{
public:
	COptimiserHennig () {}
	~COptimiserHennig () {}

	void SetParameters (const std::vector<variant>& parameters) override {};

	double GetMinimalScaleValue (std::shared_ptr<CLinearFunction> vf, CMinimiserGolden& minimiser, const arma::vec& gradient);

	arma::mat UpdateHessian () override;

	void Optimise (double target) override;
	void Optimise2 ();

private:
	arma::mat m_mVar;
	double m_rNoise;
};

//////////////////////////////////////////////////
