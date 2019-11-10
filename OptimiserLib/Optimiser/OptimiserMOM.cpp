// OptimiserPenalty.cpp
//
// 2017
//////////////////////////////////////////////////

#include "stdafx.h"

//////////////////////////////////////////////////

#include "OptimiserMOM.h"
#include "optimisable.h"
#include "CompoundMOM.h"

#include <cassert>
#include <iostream>
#include <armadillo>

#include <boost/foreach.hpp>

//////////////////////////////////////////////////

void COptimiserMOM::AddConstraint ( std::shared_ptr<const IConstraint> ptrConstraint )
{
	m_vConstraints.push_back(ptrConstraint);
	m_Multipliers.push_back (0);
}

//////////////////////////////////////////////////

void COptimiserMOM::SetParameters (const std::vector<variant>& parameters)
{
	m_rPenalty = boost::get<double> (parameters[eMOMParameters::ePenalty]);
	assert (m_rPenalty >= realZero);
}

//////////////////////////////////////////////////

void COptimiserMOM::Optimise (double target)
{
	assert (m_ptrOptimisableFunction);
	assert (m_ptrUnconstrainedOptimiser);

	// create new optimisation function using the penalty term to add the constraints
	auto multipliers = std::vector<double>();
	multipliers.resize(m_vConstraints.size()); // VS initialised the vector to zero
	auto ptrPenalisedObjective = std::make_shared<CCompoundMOMFunction>(m_ptrOptimisableFunction, m_vConstraints, multipliers);
	auto deltaOOF(realEmpty);
	double valueMainFunction(realEmpty);
	ptrPenalisedObjective->SetPenaltyFactor(m_rPenalty);

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
		ptrPenalisedObjective->SetMultipliers(m_MinVector);
		valueMainFunction = m_ptrOptimisableFunction->Evaluate (m_MinVector);// for debugging only
		++iIteration;
	} while (std::abs(deltaOOF) > m_rTolerance && iIteration <= m_uMaxIterations);

	return;
}

//////////////////////////////////////////////////