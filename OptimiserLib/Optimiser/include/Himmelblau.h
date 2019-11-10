// Himmelblau.h
//
// Typical example problem for optimiser theory
//
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "optimisable.h"

//////////////////////////////////////////////////

class FUNCTION CHimmelblau : public IFunction
{
private:
	double m_rA;
	double m_rB;
	const int m_iDimension = 2;

public:

	CHimmelblau (double rA, double rB);
	//	~CHimmelblau () {}

	double Evaluate (const arma::vec& vValues) const override;

	arma::vec CalculateGradient (const arma::vec& vec) const override;
	arma::mat CalculateHessian (const arma::vec& vec) const override;

	arma::uword GetFunctionDimension () const;
};

//////////////////////////////////////////////////

class CHimmelblauConstraint1 : public IFunction
{
private:
	const arma::uword m_iDimension = 2;

public:

	CHimmelblauConstraint1 ();
	~CHimmelblauConstraint1 () {}

	double Evaluate (const arma::vec& vValues) const override;

	arma::vec CalculateGradient (const arma::vec& vValues) const override;

	arma::mat CalculateHessian (const arma::vec& vValues) const override;

	arma::uword GetFunctionDimension () const;
};

//////////////////////////////////////////////////

class CHimmelblauConstraint2 : public IFunction // IOptimisableConstraint?
{
private:
	const arma::uword m_iDimension = 2; // needs to be set by parent 

public:

	CHimmelblauConstraint2 ();
	~CHimmelblauConstraint2 () {}

	double Evaluate (const arma::vec& vValues) const override;

	arma::vec CalculateGradient (const arma::vec& vValues) const override;

	arma::mat CalculateHessian (const arma::vec& vValues) const override;

	arma::uword GetFunctionDimension () const;
};

//////////////////////////////////////////////////