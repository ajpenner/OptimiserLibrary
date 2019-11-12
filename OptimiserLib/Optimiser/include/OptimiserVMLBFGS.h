// OptimiserVMLBFGS.h
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

namespace eLBFGSParameters
{
	enum type
	{
		eCount = 0
	};
}

namespace eLBFGSValues
{
	enum type
	{
		eDeltaGrad = 0,
		eDeltaX,
	};
};

class OPTIMISER COptimiserVMLBFGS : public BOptimiserVariableMetric
{
public:
	COptimiserVMLBFGS () {}
	~COptimiserVMLBFGS () {}

	void SetParameters (const std::vector<variant>& parameters) override;

	void Optimise (double target) override;

	arma::mat UpdateHessian () override;

	double GetMinimalScaleValue (std::shared_ptr<CLinearFunction> vf, CMinimiserGolden& minimiser, const arma::vec& gradient);

private:
	uint32_t m_Count;
	std::vector<arma::vec> m_vDeltaX;
	std::vector<arma::vec> m_vDeltaGrad;
};

//////////////////////////////////////////////////
