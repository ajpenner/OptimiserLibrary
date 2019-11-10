// Parabolic.h
//
// Simple example problem for optimiser theory
//
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "optimisable.h"

//////////////////////////////////////////////////

class FUNCTION CParabolic : public IFunction
{
private:
	double m_rA;
	double m_rB;
	double m_rC;
	const int m_iDimension = 1;

public:

	CParabolic (double rA, double rB, double rC);
	~CParabolic () {}

	double Evaluate (const arma::vec& vValues) const override;

	arma::vec CalculateGradient (const arma::vec& vec) const override;
	arma::mat CalculateHessian (const arma::vec& vec) const override;

	arma::uword GetFunctionDimension () const;
};

//////////////////////////////////////////////////