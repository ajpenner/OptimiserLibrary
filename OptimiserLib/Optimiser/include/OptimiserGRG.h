// OptimiserGRG.h
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

class OPTIMISER COptimiserGRG : public BConstraintOptimiser
{
public:
	COptimiserGRG () : m_pUnconstrainedOptimiser (nullptr), m_uMaxIterations (10), m_iSlackCount(0) {}
	~COptimiserGRG () {}

	void SetParameters (const std::vector<variant>& /*parameters*/) override {};

	void SetUnconstrainedOptimiser (IOptimiser* ptrOptimisableFunction);

	void SetInitialRadius (double rRadius) { m_rRadius = rRadius; }
	void Optimise (double target) override;

	void AddConstraint (std::shared_ptr<const IConstraint> ptrConstraint);

	void SetMaxIterations (size_t uMax) { m_uMaxIterations = uMax; }

	arma::vec FindFeasible() const;

	void SetBounds (std::vector<arma::vec>& vBounds);

private:

	void _CheckForInequalityConstraintsAndAddSlackVariables();

	IOptimiser* m_pUnconstrainedOptimiser;

	arma::vec m_vCoordinateMin;
	arma::vec m_vCoordinateMax;

	arma::vec m_mDerivatives;

	double m_rRadius;
	std::vector<std::shared_ptr<const IConstraint>> m_vConstraints;

	std::vector<arma::vec> m_vCoordinateBounds;

	size_t m_uMaxIterations;

	arma::uword m_iSlackCount;
};

//////////////////////////////////////////////////