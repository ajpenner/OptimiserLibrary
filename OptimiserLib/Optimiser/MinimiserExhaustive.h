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

class OPTIMISER CMinimiserExhaustive : public BExtremum
{
public:
	CMinimiserExhaustive () {}
	~CMinimiserExhaustive () {}

	void SetParameters (const std::vector<variant>& parameters) override {};

	void Optimise (double /*target*/) override;

private:

	void CheckAndSwap (double& rLeft, double& rRight);

	int m_iIntermediatePoints;
};