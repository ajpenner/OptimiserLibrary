// OptimiserNewton.h
//
// 2016
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "optimiser.h"
#include "optimisable.h"
#include <limits>
#include <armadillo>

//////////////////////////////////////////////////

class OPTIMISER COptimiserNewton : public BOptimiser
{
public:
	COptimiserNewton () {}
	~COptimiserNewton () {}

	void SetLeftBracket (const arma::vec& /*vec*/) override { return; }
	void SetRightBracket (const arma::vec& /*vec*/) override { return; }

	void SetParameters (const std::vector<variant>& /*parameters*/) override {};

	void Optimise (double target) override;

private:

};

//////////////////////////////////////////////////