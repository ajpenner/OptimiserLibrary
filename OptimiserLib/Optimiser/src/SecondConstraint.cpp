// FirstConstraint.cpp
//
// 2017
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif
#include "SecondConstraint.h"
#include <cassert>

//////////////////////////////////////////////////

CSecondConstraint::CSecondConstraint () : m_rA(std::numeric_limits<double>::max()), m_rB(std::numeric_limits<double>::max())
{
}

//////////////////////////////////////////////////

void CSecondConstraint::SetParameters(const arma::vec& vector)
{
	assert (vector.size () >= 2);
	m_rA = vector.at (0);
	m_rB = vector.at (1);
	m_iDimension = 2;
}

//////////////////////////////////////////////////

double CSecondConstraint::Evaluate(const arma::vec& vValues) const
{
	assert (vValues.size() == 2);
	return m_rA + m_rB*vValues.at(0) - vValues.at(1);
}

//////////////////////////////////////////////////

arma::vec CSecondConstraint::CalculateGradient (const arma::vec& vec) const
{ 
	return arma::vec (m_iDimension).zeros (); 
}

//////////////////////////////////////////////////

arma::mat CSecondConstraint::CalculateHessian (const arma::vec& vec) const
{ 
	return arma::mat (m_iDimension, m_iDimension).zeros (); 
}

//////////////////////////////////////////////////

arma::uword CSecondConstraint::GetFunctionDimension () const
{ 
	return m_iDimension; 
}

//////////////////////////////////////////////////
