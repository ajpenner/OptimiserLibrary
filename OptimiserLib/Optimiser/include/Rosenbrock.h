// Rosenbrock.h
//
// Typical example problem for optimiser theory
//
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "optimisable.h"

//////////////////////////////////////////////////

class FUNCTION CRosenbrock : public IFunction
{
private:
	double m_rA;
	double m_rB;
	const int m_iDimension = 2;

public:

	CRosenbrock (double rA, double rB);
	~CRosenbrock () {}

	double Evaluate (const arma::vec& vValues) const override;

	arma::vec CalculateGradient (const arma::vec& vec) const override;

	arma::mat CalculateHessian (const arma::vec& vec) const override;

	arma::uword GetFunctionDimension () const;
};

//////////////////////////////////////////////////