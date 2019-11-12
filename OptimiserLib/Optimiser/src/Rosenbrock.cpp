// Rosenbrock.cpp
//
// Typical example problem for optimiser theory
//
// 2017
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif
#include "Rosenbrock.h"
#include <cassert>

//////////////////////////////////////////////////

CRosenbrock::CRosenbrock (double rA, double rB) : m_rA (rA), m_rB (rB)
{
}

//////////////////////////////////////////////////

double CRosenbrock::Evaluate (const arma::vec& vValues) const
{
	assert (vValues.size () == m_iDimension);
	auto rX = vValues (0);
	auto rY = vValues (1);
	return std::pow ((m_rA - rX), 2) + m_rB*std::pow ((rY - rX*rX), 2);
}

//////////////////////////////////////////////////

arma::vec CRosenbrock::CalculateGradient (const arma::vec& vec) const
{
	assert (vec.size () == m_iDimension);
	auto vGrad = arma::vec (m_iDimension);
	vGrad (0) = -2.0*(m_rA - vec (0)) - 4.0*m_rB*vec (0)*(vec (1) - vec (0)*vec (0));
	vGrad (1) = 2.0*m_rB*(vec (1) - vec (0) * vec (0));
	return vGrad;
};

//////////////////////////////////////////////////

arma::mat CRosenbrock::CalculateHessian (const arma::vec& vec) const
{
	assert (vec.size () == m_iDimension);
	auto mHess = arma::mat (m_iDimension, m_iDimension);
	mHess (0, 0) = 2.0*m_rA - 4.0*m_rB*(vec (1) - vec (0) * vec (0)) + 8.0*m_rB*vec (0) * vec (0);
	mHess (0, 1) = -4.0*m_rB * vec (0);
	mHess (1, 0) = -4.0*m_rB * vec (0);
	mHess (1, 1) = 2.0*m_rB;
	return mHess;
}

//////////////////////////////////////////////////

arma::uword CRosenbrock::GetFunctionDimension () const
{
	return m_iDimension;
}

//////////////////////////////////////////////////
