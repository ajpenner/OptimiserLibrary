// OptimiserVMBroyden.h
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

namespace eBroydenParameters
{
	enum type
	{
		ePhi = 0
	};
}

class OPTIMISER COptimiserVMBroyden : public BOptimiserVariableMetric
{
public:
	COptimiserVMBroyden() : m_phi(0.5) {}
	~COptimiserVMBroyden() {}

	void SetParameters (const std::vector<variant>& parameters) override;

	arma::mat UpdateHessian () override;

	void Optimise(double target) override;

	double GetMinimalScaleValue (std::shared_ptr<CLinearFunction> vf, CMinimiserGolden& minimiser, const arma::vec& gradient);

private:
	double m_phi;
};

//////////////////////////////////////////////////