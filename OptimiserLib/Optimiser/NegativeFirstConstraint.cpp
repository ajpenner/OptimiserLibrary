// FirstConstraint.cpp
//
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "stdafx.h"
#include "NegativeFirstConstraint.h"

//////////////////////////////////////////////////

CNegativeFirstConstraint::CNegativeFirstConstraint () : m_rA(realEmpty), m_rB(realEmpty)
{
}

//////////////////////////////////////////////////

void CNegativeFirstConstraint::SetParameters(const arma::vec& vector)
{
	assert (vector.size () >= 2);
	m_rA = vector.at (0);
	m_rB = vector.at (1);
	m_iDimension = 2;
}

//////////////////////////////////////////////////

double CNegativeFirstConstraint::Evaluate(const arma::vec& vValues) const
{
	assert (vValues.size() == 2);
	return -std::pow ((vValues.at (0) + m_rA), 2) - std::pow (vValues.at (1), 2) + m_rB;
}

//////////////////////////////////////////////////

arma::vec CNegativeFirstConstraint::CalculateGradient (const arma::vec& vec) const
{ 
	return arma::vec (m_iDimension).zeros (); 
}

//////////////////////////////////////////////////

arma::mat CNegativeFirstConstraint::CalculateHessian (const arma::vec& vec) const
{ 
	return arma::mat (m_iDimension, m_iDimension).zeros (); 
}

//////////////////////////////////////////////////

arma::uword CNegativeFirstConstraint::GetFunctionDimension () const
{ 
	return m_iDimension; 
}

//////////////////////////////////////////////////