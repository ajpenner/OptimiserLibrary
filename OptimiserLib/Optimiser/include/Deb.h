// Deb.h
//
// Simple example problem for optimiser theory
//
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "optimisable.h"

//////////////////////////////////////////////////

class FUNCTION CDeb : public IFunction
{
private:
	double m_rA;
	const int m_iDimension = 1;

public:

	CDeb (double rA);
	~CDeb () {}

	double Evaluate (const arma::vec& vValues) const override;

	arma::vec CalculateGradient (const arma::vec& vec) const override;
	arma::mat CalculateHessian (const arma::vec& vec) const override;

	arma::uword GetFunctionDimension () const;
};

//////////////////////////////////////////////////