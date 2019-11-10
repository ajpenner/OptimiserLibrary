// MinimiserBrent.h
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

class OPTIMISER CMinimiserBrent : public BExtremum
{
public:
	CMinimiserBrent () {}
	~CMinimiserBrent () {}

	void SetParameters (const std::vector<variant>& parameters) override {};

	void Optimise (double /*target*/) override;

private:

	void _CheckAndSwap (double& rLeft, double& rRight);
	void _PerformSearch (double& rLeft, double& rRight);

};

//////////////////////////////////////////////////