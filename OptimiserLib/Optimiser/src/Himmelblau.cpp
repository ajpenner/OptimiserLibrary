// Himmelblau.h
//
// Typical example problem for optimiser theory
//
// 2017
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif
#include "Himmelblau.h"
#include <cassert>

//////////////////////////////////////////////////

CHimmelblau::CHimmelblau (double rA, double rB) : m_rA (rA), m_rB (rB)
{
}

//////////////////////////////////////////////////

double CHimmelblau::Evaluate (const arma::vec& vValues) const
{
	assert (vValues.size () == m_iDimension);
	auto rX1 = vValues (0);
	auto rX2 = vValues (1);
	auto rPart1 = std::pow (rX1, 2) + rX2 + m_rA;
	auto rPart2 = std::pow (rX2, 2) + rX1 + m_rB;

	return std::pow(rPart1,2) + std::pow(rPart2,2);
}

//////////////////////////////////////////////////

arma::vec CHimmelblau::CalculateGradient (const arma::vec& vec) const
{
	assert (vec.size () == m_iDimension);
	auto& x1 = vec.at(0);
	auto& x2 = vec.at(1);
	auto rPart1 = std::pow (x1, 2) + x2 + m_rA;
	auto rPart2 = std::pow (x2, 2) + x1 + m_rB;

	double f1 = 4 * x1 * rPart1 + 2 * rPart2;
	double f2 = 2 * rPart1 + 4 * x2 * rPart2;
	return arma::vec({f1,f2});
}

//////////////////////////////////////////////////

arma::mat CHimmelblau::CalculateHessian (const arma::vec& vec) const
{ 
	assert (vec.size () == m_iDimension);
	auto mHess = arma::mat (m_iDimension, m_iDimension).zeros ();
	auto& x1 = vec.at (0);
	auto& x2 = vec.at (1);
	auto rPart1 = std::pow (x1, 2) + x2 + m_rA;
	auto rPart2 = std::pow (x2, 2) + x1 + m_rB;

	mHess(0, 0) = 2 + 8*std::pow(x1,2) + 4*(rPart1);
	mHess(1, 0) = 4 * x1 + 4 * x2;
	mHess(0, 1) = mHess(1, 0);
	mHess(1, 1) = 2 + 8 * std::pow (x2, 2) + 4 * (rPart2);
	return mHess;
}

//////////////////////////////////////////////////

arma::uword CHimmelblau::GetFunctionDimension () const 
{ 
	return m_iDimension; 
}

//////////////////////////////////////////////////

CHimmelblauConstraint1::CHimmelblauConstraint1 ()
{
}

//////////////////////////////////////////////////

double CHimmelblauConstraint1::Evaluate (const arma::vec& vValues) const
{
	assert (vValues.size () >= m_iDimension);
	return 26.0 - std::pow (vValues[0] - 5.0, 2) - std::pow (vValues[1], 2);
}

//////////////////////////////////////////////////

arma::vec CHimmelblauConstraint1::CalculateGradient (const arma::vec& vValues) const
{
	auto vReturn = arma::vec (m_iDimension).zeros ();
	vReturn[0] = -2.0*(vValues[0] - 5.0);
	vReturn[1] = -2.0*vValues[1];
	return vReturn;
}

//////////////////////////////////////////////////

arma::mat CHimmelblauConstraint1::CalculateHessian (const arma::vec& vValues) const
{
	auto mReturn = arma::mat (m_iDimension, m_iDimension).zeros ();
	mReturn (0, 0) = -2.0;
	mReturn (1, 1) = -2.0;
	return mReturn;
}

//////////////////////////////////////////////////

arma::uword CHimmelblauConstraint1::GetFunctionDimension () const
{
	return m_iDimension;
}

//////////////////////////////////////////////////

CHimmelblauConstraint2::CHimmelblauConstraint2 ()
{
}

//////////////////////////////////////////////////

double CHimmelblauConstraint2::Evaluate (const arma::vec& vValues) const
{
	assert (vValues.size () >= m_iDimension);
	return 20.0 - 4.0*vValues[0] - vValues[1];
}

//////////////////////////////////////////////////

arma::vec CHimmelblauConstraint2::CalculateGradient (const arma::vec& vValues) const
{
	auto vReturn = arma::vec (m_iDimension).zeros ();
	vReturn[0] = -4.0;
	vReturn[1] = -1.0;
	return vReturn;
}

//////////////////////////////////////////////////

arma::mat CHimmelblauConstraint2::CalculateHessian (const arma::vec& vValues) const
{
	auto mReturn = arma::mat (m_iDimension, m_iDimension).zeros ();
	mReturn (0, 0) = 0x0;
	mReturn (1, 1) = 0x0;
	return mReturn;
}

//////////////////////////////////////////////////

arma::uword CHimmelblauConstraint2::GetFunctionDimension () const
{
	return m_iDimension;
}

//////////////////////////////////////////////////
