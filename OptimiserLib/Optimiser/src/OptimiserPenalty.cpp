// OptimiserPenalty.cpp
//
// 2017
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif

//////////////////////////////////////////////////

#include "OptimiserPenalty.h"
#include "optimisable.h"
#include "CompoundPenalty.h"

#include <cassert>
#include <armadillo>

#include <boost/foreach.hpp>

//////////////////////////////////////////////////

void COptimiserPenalty::AddConstraint ( std::shared_ptr<const IConstraint> ptrConstraint )
{
	m_vConstraints.push_back (ptrConstraint);
}

//////////////////////////////////////////////////

void COptimiserPenalty::SetParameters (const std::vector<variant>& parameters)
{
	m_rPenalty = boost::get<double> (parameters[ePenaltyParameters::ePenalty]);
	assert (m_rPenalty >= 0x0);
	m_rScaleFactor = boost::get<double> (parameters[ePenaltyParameters::eScaleFactor]);
	assert (m_rScaleFactor >= 1);
}

//////////////////////////////////////////////////

void COptimiserPenalty::Optimise (double target)
{
	assert (m_ptrOptimisableFunction);
	assert (m_ptrUnconstrainedOptimiser);

	// create new optimisation function using the penalty term to add the constraints
	auto coeff = std::vector<double>();
	coeff.resize(m_vConstraints.size());
	std::generate_n(coeff.begin(), m_vConstraints.size(), [&]() { return m_rPenalty; });
	auto ptrPenalisedObjective = std::make_shared<CCompoundPenaltyFunction>(m_ptrOptimisableFunction, m_vConstraints, coeff);
	auto deltaOOF(std::numeric_limits<double>::max());
	double valueMainFunction(std::numeric_limits<double>::max());
	size_t iIteration (0);
	do
	{
		deltaOOF = std::abs(m_ptrUnconstrainedOptimiser->GetOOF());
		m_ptrUnconstrainedOptimiser->SetOptimisableFunction(ptrPenalisedObjective);
		m_ptrUnconstrainedOptimiser->SetInitialGuess(m_MinVector);
		m_ptrUnconstrainedOptimiser->SetTolerance(m_rTolerance);
		m_ptrUnconstrainedOptimiser->Optimise();
		m_MinVector = m_ptrUnconstrainedOptimiser->GetVector();
		m_OOFValue = m_ptrUnconstrainedOptimiser->GetOOF();
		deltaOOF -= std::abs(m_OOFValue);
		m_rPenalty *= m_rScaleFactor;
		std::transform (coeff.begin(), coeff.end(), coeff.begin(), [&](auto& x) { return m_rPenalty; });
		ptrPenalisedObjective->SetCoeff (coeff);
		valueMainFunction = m_ptrOptimisableFunction->Evaluate (m_MinVector);
		++iIteration;
	} while (std::abs(deltaOOF) > m_rTolerance && iIteration <= m_uMaxIterations);

	return;
}

//////////////////////////////////////////////////
