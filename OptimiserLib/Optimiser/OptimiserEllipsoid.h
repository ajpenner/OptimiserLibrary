// OptimiserEllipsoid.h
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

namespace eEllipsoidParameters
{
	enum type
	{
		eRadius = 0
	};
}

class OPTIMISER COptimiserEllipsoid : public BConstraintOptimiser
{
public:
	COptimiserEllipsoid () : m_pUnconstrainedOptimiser(nullptr), m_uMaxIterations(10) {}
	~COptimiserEllipsoid () {}

	void SetParameters (const std::vector<variant>& parameters) override;

	void SetUnconstrainedOptimiser (IOptimiser* ptrOptimisableFunction);

	void Optimise (double target) override;

	void AddConstraint (std::shared_ptr<const IConstraint> ptrConstraint) override;

private:

	IOptimiser* m_pUnconstrainedOptimiser;

	double m_rRadius;
	std::vector<std::shared_ptr<const ILinearConstraint>> m_vConstraints; // Just IConstraint, no reason to be restrictive

	size_t m_uMaxIterations;
};

//////////////////////////////////////////////////