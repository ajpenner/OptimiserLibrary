// COptimiserHJ.cpp
// Hooke - Jeeves
// 2017
//////////////////////////////////////////////////

#include "stdafx.h"

//////////////////////////////////////////////////

#include "OptimiserHJ.h"
#include "optimisable.h"
#include <cassert>
#include <iostream>

//////////////////////////////////////////////////

COptimiserHJ::COptimiserHJ () :
	m_IncrementFactor(0.5),
	m_ReductionFactor(0.5)
{ 
}

//////////////////////////////////////////////////

void COptimiserHJ::SetLeftBracket (const arma::vec& vLeft)
{
	m_LeftBracket = vLeft; 
}

//////////////////////////////////////////////////

void COptimiserHJ::SetRightBracket (const arma::vec& vRight)
{ 
	m_RightBracket = vRight; 
}

//////////////////////////////////////////////////

double COptimiserHJ::GetInverseOOF() const
{
	return m_InverseOOFValue;
}

//////////////////////////////////////////////////

void COptimiserHJ::InverseOptimise(double target)
{
	Optimise(-target);
}

//////////////////////////////////////////////////

void COptimiserHJ::SetParameters (const std::vector<variant>& parameters)
{
	assert (parameters.size () >= 2);
	m_IncrementFactor = boost::get<double>(parameters[eHJParameters::eIncrement]);
	m_ReductionFactor = boost::get<double>(parameters[eHJParameters::eReduction]);
};

//////////////////////////////////////////////////

void COptimiserHJ::CreateSimplex() // I do not understand this yet
{

}

//////////////////////////////////////////////////

arma::vec COptimiserHJ::Explore(arma::vec& start, bool& bSuccess)
{
	auto eval = std::vector<double> (m_ptrOptimisableFunction->GetFunctionDimension () + 1);
	auto coord = std::vector<arma::vec> (m_ptrOptimisableFunction->GetFunctionDimension () + 1);

	// explore
	auto cStart = start;
	do
	{
		for (auto index = 0; index < m_ptrOptimisableFunction->GetFunctionDimension (); ++index)
		{
			eval[0] = m_ptrOptimisableFunction->Evaluate (cStart);
			coord[0] = cStart; // if coord[0] points to the memory location of m_MinVector we can skip this assignment

			auto delta = m_I.col (index)*m_IncrementFactor;
			coord[1] = cStart + delta;
			eval[1] = m_ptrOptimisableFunction->Evaluate (cStart + delta);
			coord[2] = cStart - delta;
			eval[2] = m_ptrOptimisableFunction->Evaluate (cStart - delta);
			auto itEval = std::min_element (eval.cbegin (), eval.cend ());
			cStart = coord[std::distance (eval.cbegin (), itEval)];
		}
		if (!arma::norm(cStart - start))
		{
			m_IncrementFactor *= m_ReductionFactor;
		}
	} while(!arma::norm(cStart - start) && m_IncrementFactor > m_rTolerance);

	if (m_IncrementFactor < m_rTolerance)
		bSuccess = false;

	return cStart;
}

//////////////////////////////////////////////////

void COptimiserHJ::Optimise (double target) // not completed
{
	assert (m_ptrOptimisableFunction);
	m_I = arma::eye (m_ptrOptimisableFunction->GetFunctionDimension (), m_ptrOptimisableFunction->GetFunctionDimension ());

	m_OOFValue = realEmpty;
	arma::vec x2 = arma::zeros(m_ptrOptimisableFunction->GetFunctionDimension());
	arma::vec xp = m_MinVector;
	bool bSuccess = true;
	do
	{
		x2 = Explore(xp, bSuccess);
		xp = (x2 + (x2 - m_MinVector));

		m_MinVector = x2;
		x2 = Explore(xp, bSuccess);
		auto eval1 = m_ptrOptimisableFunction->Evaluate(x2);
		auto eval2 = m_ptrOptimisableFunction->Evaluate(m_MinVector); // make this m_OOFValue
		if (eval1 < eval2)
		{
			xp = (x2 + (x2 - m_MinVector));
		}
		else if(m_IncrementFactor > m_rTolerance)
		{
			m_IncrementFactor *= m_ReductionFactor;
			x2 = Explore(m_MinVector, bSuccess);
		}

		if(!bSuccess)
		{
			m_OOFValue = m_ptrOptimisableFunction->Evaluate(x2);
		}
		
		m_MinVector = x2;
	} while (m_OOFValue - target > m_rTolerance);
		
	return;
}

//////////////////////////////////////////////////