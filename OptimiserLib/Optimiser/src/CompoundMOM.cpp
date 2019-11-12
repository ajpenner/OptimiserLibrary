// Compound.cpp
//
// Simple example problem for optimiser theory
//
// 2017
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif
#include "CompoundMOM.h"
#include <memory>
#include "boost/iterator/zip_iterator.hpp"

//////////////////////////////////////////////////

CCompoundMOMFunction::CCompoundMOMFunction (std::shared_ptr<const IFunction> ptrOOF, std::vector<std::shared_ptr<const IConstraint>> vConstraints, std::vector<double> vMultipliers) :
	m_ptrOOF (ptrOOF), m_Constraints (vConstraints), m_Multiplier(vMultipliers), m_rPenalty(0x0)
{
	assert (m_Multiplier.size () == m_Constraints.size ());
}

//////////////////////////////////////////////////

double CCompoundMOMFunction::Evaluate(const arma::vec& vValues) const
{
	double value = m_ptrOOF->Evaluate (vValues);

	auto begin = boost::make_zip_iterator (	boost::make_tuple (m_Constraints.begin(), m_Multiplier.begin()) );
	auto end = boost::make_zip_iterator (boost::make_tuple (m_Constraints.end (), m_Multiplier.end ()));

	double rPenalty = m_rPenalty; // should be able to pass a member variable to a lambda...
	std::for_each (begin, end, 
		[&value, rPenalty, vValues](const boost::tuple<const std::shared_ptr<const IConstraint>&, const double&>& t)
	{
		auto constraint = t.get<eMOMCompoundParameters::eConstraint> ();
		auto multiplier = t.get<eMOMCompoundParameters::eMultiplier> ();

		auto constraintViolation = constraint->Evaluate (vValues) + multiplier;
		if (constraintViolation > 0) constraintViolation = 0;
		
		value += rPenalty*(std::pow (constraintViolation, 2) - std::pow (multiplier, 2));
	});

	return value;
}

//////////////////////////////////////////////////

arma::vec CCompoundMOMFunction::CalculateGradient (const arma::vec& vec) const
{
	arma::vec value = m_ptrOOF->CalculateGradient (vec);
	for(const std::shared_ptr<const IConstraint>& var : m_Constraints)
	{
		value += var->CalculateGradient (vec);
	}
	return value;
}

//////////////////////////////////////////////////

arma::mat CCompoundMOMFunction::CalculateHessian (const arma::vec& vec) const
{
	arma::mat value = m_ptrOOF->CalculateHessian (vec);
	for (const std::shared_ptr<const IConstraint>& var : m_Constraints)
	{
		value += var->CalculateHessian (vec);
	}
	return value;
}

//////////////////////////////////////////////////

arma::uword CCompoundMOMFunction::GetFunctionDimension () const
{
	return m_ptrOOF->GetFunctionDimension ();
};

//////////////////////////////////////////////////

std::vector<double> CCompoundMOMFunction::GetConstraintViolations (const arma::vec& vec) const
{
	assert (vec.size () == m_ptrOOF->GetFunctionDimension ());
	auto vMultipliers = std::vector<double> ();
	std::transform (m_Constraints.cbegin(), m_Constraints.cend (), std::back_inserter(vMultipliers), 
		[vec](const auto constraint )
	{
		return constraint->Evaluate(vec);
	});

	return vMultipliers;
}

//////////////////////////////////////////////////

void CCompoundMOMFunction::SetMultipliers (const arma::vec& vec)
{
	assert (vec.size () == m_ptrOOF->GetFunctionDimension ());
	assert (m_Constraints.size () == m_Multiplier.size());

	std::transform (m_Constraints.cbegin (), m_Constraints.cend (), m_Multiplier.begin(),
		[vec](const auto constraint)
	{
		auto violation = 0.0;
		if(!constraint->IsSatisfied(vec))
			violation = constraint->Evaluate (vec);
		return violation;
	});
}

//////////////////////////////////////////////////

void CCompoundMOMFunction::SetPenaltyFactor( double penalty )
{
	m_rPenalty = penalty;
}

//////////////////////////////////////////////////
