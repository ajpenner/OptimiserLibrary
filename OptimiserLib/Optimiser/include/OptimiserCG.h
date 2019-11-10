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

class OPTIMISER COptimiserCG : public BOptimiser
{
public:
	COptimiserCG () {}
	~COptimiserCG () {}

	void SetLeftBracket (const arma::vec& /*vec*/) override { return; }
	void SetRightBracket (const arma::vec& /*vec*/) override { return; }

	void SetParameters (const std::vector<variant>& /*parameters*/) override {};

	void Optimise (double target) override;

	double GetMinimalScaleValue (std::shared_ptr<CLinearFunction> vf, CMinimiserGolden& minimiser, const arma::vec& gradient);

private:

};

//////////////////////////////////////////////////