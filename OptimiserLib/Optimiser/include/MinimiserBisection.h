// MinimiserBisection.h
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

class OPTIMISER CMinimiserBisection : public BExtremum
{
public:
	CMinimiserBisection () {}
	~CMinimiserBisection () {}

	void SetParameters (const std::vector<variant>& parameters) override {};
	void Optimise (double /*target*/) override;

private:

};

//////////////////////////////////////////////////