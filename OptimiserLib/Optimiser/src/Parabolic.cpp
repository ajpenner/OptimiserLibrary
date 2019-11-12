// Parabolic.cpp
//
// Simple example problem for optimiser theory
//
// 2017
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif
#include "Parabolic.h"
#include <cassert>

//////////////////////////////////////////////////

CParabolic::CParabolic (double rA, double rB, double rC) : m_rA (rA), m_rB (rB), m_rC (rC) 
{
}

//////////////////////////////////////////////////

double CParabolic::Evaluate(const arma::vec& vValues) const
{
	assert (vValues.n_cols == m_iDimension);
	auto rX = vValues (0);
	return m_rA*std::pow (rX, 2) + m_rB*rX + m_rC;	
}

//////////////////////////////////////////////////

arma::vec CParabolic::CalculateGradient (const arma::vec& vec) const 
{ 
	return arma::vec (m_iDimension).zeros (); 
}

//////////////////////////////////////////////////

arma::mat CParabolic::CalculateHessian (const arma::vec& vec) const 
{ 
	return arma::mat (m_iDimension, m_iDimension).zeros (); 
}

//////////////////////////////////////////////////

arma::uword CParabolic::GetFunctionDimension () const 
{ 
	return m_iDimension; 
}

//////////////////////////////////////////////////
