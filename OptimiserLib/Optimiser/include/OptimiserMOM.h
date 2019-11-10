// OptimiserMOM.h
// Method of Multipliers
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "optimiser.h"
#include "optimisable.h"
#include "constraints.h"
#include "CompoundMOM.h"

#include <limits>
#include <armadillo>

//////////////////////////////////////////////////

namespace eMOMParameters
{
	enum type
	{
		ePenalty = 0
	};
}

class OPTIMISER COptimiserMOM : public BConstraintOptimiser
{
public:
	COptimiserMOM () {}
	~COptimiserMOM () {}

	void SetParameters (const std::vector<variant>& parameters) override;

	void Optimise (double target) override;

	void AddConstraint (std::shared_ptr<const IConstraint> ptrConstraint) override;

private:

	double m_rPenalty;
	std::vector<double> m_Multipliers;
	std::vector<std::shared_ptr<const IConstraint>> m_vConstraints;
};

//////////////////////////////////////////////////