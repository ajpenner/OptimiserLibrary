// NegativeFirstConstraint.h
//
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "optimisable.h"
#include <memory>

//////////////////////////////////////////////////

class FUNCTION CNegativeFirstConstraint : public IFunction
{
private:
	uint8_t m_iDimension;
	double m_rA;
	double m_rB;

public:

	CNegativeFirstConstraint ();
	~CNegativeFirstConstraint () {}

	void SetParameters(const arma::vec& vector);

	double Evaluate (const arma::vec& vValues) const override;

	arma::vec CalculateGradient (const arma::vec& vec) const override;
	arma::mat CalculateHessian (const arma::vec& vec) const override;

	arma::uword GetFunctionDimension () const;
};

//////////////////////////////////////////////////