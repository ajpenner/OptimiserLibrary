// Compound.cpp
//
// Simple example problem for optimiser theory
//
// 2017
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif
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
	for (const std::shared_ptr<const IConstraint>& var : m_Constraints)
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
	for (const std::shared_ptr<const IConstraint>& var : m_Constraints)
	{
		value += m_Coeff[0] * var->CalculateGradient (vec);
	}
	return value;
}

//////////////////////////////////////////////////

arma::mat CCompoundPenaltyFunction::CalculateHessian (const arma::vec& vec) const
{
	arma::mat value = m_ptrOOF->CalculateHessian (vec);
	for (const std::shared_ptr<const IConstraint>& var : m_Constraints)
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
