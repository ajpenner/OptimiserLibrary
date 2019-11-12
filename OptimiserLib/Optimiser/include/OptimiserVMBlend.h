// OptimiserVMBroyden.h
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

namespace eBlendParameters
{
	enum type
	{
		eBlendRatio = 0
	};
}

class OPTIMISER COptimiserVMBlend : public BOptimiserVariableMetric
{
public:
	COptimiserVMBlend () : m_blendRatio (0.5) {}
	~COptimiserVMBlend () {}

	void SetParameters (const std::vector<variant>& parameters) override;

	arma::mat UpdateHessian () override;

	void Optimise(double target) override;

	double GetMinimalScaleValue (std::shared_ptr<CLinearFunction> vf, CMinimiserGolden& minimiser, const arma::vec& gradient);

private:
	double m_blendRatio;
	arma::mat UpdateHessianBFGS();
	arma::mat UpdateHessianDFP();

};

//////////////////////////////////////////////////
