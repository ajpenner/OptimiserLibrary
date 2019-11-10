// MinimiserBoxed.h
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

class OPTIMISER CMinimiserBoxed : public BExtremum
{
public:
	CMinimiserBoxed () {}
	CMinimiserBoxed ( BExtremum* pcUnconstrainedOptimiser ) {}
	~CMinimiserBoxed () {}

	void SetParameters (const std::vector<variant>& parameters) override {};

	void Optimise (double /*target*/) override;

private:

};

//////////////////////////////////////////////////