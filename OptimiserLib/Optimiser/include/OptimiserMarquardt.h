// OptimiserMarquardt.h
//
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "optimiser.h"
#include "optimisable.h"

#include "LinearFunction.h"
#include "MinimiserGolden.h"

#include <limits>
#include <armadillo>

//////////////////////////////////////////////////

namespace eMarquardtParameters
{
	enum type
	{
		ePenalty = 0
	};
}

class OPTIMISER COptimiserMarquardt : public BOptimiser
{
public:
	COptimiserMarquardt () : m_penalty(10000) {}
	~COptimiserMarquardt () {}

	void SetParameters (const std::vector<variant>& parameters) override;

	void Optimise (double target) override;

	double GetMinimalScaleValue (std::shared_ptr<CLinearFunction> ptrVF, CMinimiserGolden& minimiser, const arma::vec& gradient);

private:

	double m_penalty;
};

//////////////////////////////////////////////////