// Deb.cpp
//
// Simple example problem for optimiser theory
//
// 2017
//////////////////////////////////////////////////

#ifdef _WIN32 
#include "stdafx.h"
#endif
#include "Deb.h"
#include <cassert>

//////////////////////////////////////////////////

CDeb::CDeb (double rA) : m_rA (rA)
{
}

//////////////////////////////////////////////////

double CDeb::Evaluate(const arma::vec& vValues) const
{
	assert (vValues.n_cols == m_iDimension);
	auto rX = vValues (0);
	return std::pow (rX, 2) + m_rA/rX;	
}

//////////////////////////////////////////////////

arma::vec CDeb::CalculateGradient (const arma::vec& vec) const
{ 
	return arma::vec(); 
}

//////////////////////////////////////////////////

arma::mat CDeb::CalculateHessian (const arma::vec& vec) const
{ 
	return arma::mat(); 
}

//////////////////////////////////////////////////

arma::uword CDeb::GetFunctionDimension () const
{ 
	return m_iDimension; 
}

//////////////////////////////////////////////////
