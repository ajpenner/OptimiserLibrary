// FirstConstraint.cpp
//
// 2017
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif
#include "FirstConstraint.h"
#include <cassert>

//////////////////////////////////////////////////

CFirstConstraint::CFirstConstraint () : m_rA(std::numeric_limits<double>::max()), m_rB(std::numeric_limits<double>::max())
{
}

//////////////////////////////////////////////////

void CFirstConstraint::SetParameters(const arma::vec& vector)
{
	assert (vector.size () >= 2);
	m_rA = vector.at (0);
	m_rB = vector.at (1);
	m_iDimension = 2;
}

//////////////////////////////////////////////////

double CFirstConstraint::Evaluate(const arma::vec& vValues) const
{
	assert (vValues.size() == 2);
	return std::pow ((vValues.at (0) + m_rA), 2) + std::pow (vValues.at (1), 2) + m_rB;
}

//////////////////////////////////////////////////

arma::vec CFirstConstraint::CalculateGradient (const arma::vec& vec) const
{ 
	return arma::vec (m_iDimension).zeros (); 
}

//////////////////////////////////////////////////

arma::mat CFirstConstraint::CalculateHessian (const arma::vec& vec) const
{ 
	return arma::mat (m_iDimension, m_iDimension).zeros (); 
}

//////////////////////////////////////////////////

arma::uword CFirstConstraint::GetFunctionDimension () const
{ 
	return m_iDimension; 
}

//////////////////////////////////////////////////
