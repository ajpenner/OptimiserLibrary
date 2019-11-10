// OptimiserPenalty.h
//
// 2016
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "optimiser.h"
#include "optimisable.h"
#include "constraints.h"

#include <limits>
#include <armadillo>

//////////////////////////////////////////////////

namespace ePenaltyParameters
{
	enum type
	{
		ePenalty = 0,
		eScaleFactor
	};
}

class OPTIMISER COptimiserPenalty : public BConstraintOptimiser
{
public:
	COptimiserPenalty () {}
	~COptimiserPenalty () {}

	void SetParameters (const std::vector<variant>& parameters) override;

	void Optimise (double target) override;

	void AddConstraint (std::shared_ptr<const IConstraint> ptrConstraint) override;

private:

	double m_rPenalty;
	double m_rScaleFactor;
	std::vector<std::shared_ptr<const IConstraint>> m_vConstraints;
};

//////////////////////////////////////////////////