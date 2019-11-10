// Compound.cpp
//
// Simple example problem for optimiser theory
//
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "stdafx.h"
#include "CompoundPenalty.h"
#include <memory>

//////////////////////////////////////////////////

CCompoundPenaltyFunction::CCompoundPenaltyFunction (std::shared_ptr<const IFunction> ptrOOF, std::vector<std::shared_ptr<const IConstraint>> vConstraints, std::vector<double> vCoeff) :
	m_ptrOOF (ptrOOF), m_Constraints (vConstraints), m_Coeff (vCoeff)
{
	assert (m_Coeff.size () == m_Constraints.size ());
}

//////////////////////////////////////////////////

double CCompoundPenaltyFunction::Evaluate(const arma::vec& vValues) const
{
	double value = m_ptrOOF->Evaluate (vValues);
	for each (const std::shared_ptr<const IConstraint> var in m_Constraints)
	{
		auto rPenalty = m_Coeff[0] * var->Evaluate (vValues);
		if (!var->IsSatisfied (vValues))
			value += std::pow (rPenalty, 2);
	}
	return value;
}

//////////////////////////////////////////////////

arma::vec CCompoundPenaltyFunction::CalculateGradient (const arma::vec& vec) const
{
	arma::vec value = m_ptrOOF->CalculateGradient (vec);
	for each (const std::shared_ptr<const IConstraint> var in m_Constraints)
	{
		value += m_Coeff[0] * var->CalculateGradient (vec);
	}
	return value;
}

//////////////////////////////////////////////////

arma::mat CCompoundPenaltyFunction::CalculateHessian (const arma::vec& vec) const
{
	arma::mat value = m_ptrOOF->CalculateHessian (vec);
	for each (const std::shared_ptr<const IConstraint> var in m_Constraints)
	{
		value += m_Coeff[0] * var->CalculateHessian (vec);
	}
	return value;
}

//////////////////////////////////////////////////

arma::uword CCompoundPenaltyFunction::GetFunctionDimension () const
{
	return m_ptrOOF->GetFunctionDimension ();
};

//////////////////////////////////////////////////