// constraints.h
//
// 2016
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "Optimisable.h"
#include <armadillo> // should forward declare

//////////////////////////////////////////////////

enum EConstraintType
{
	eLT,
	eLE,
	eGT,
	eGE,
	eEQ,

	eUnknown,
};

//////////////////////////////////////////////////

class FUNCTION IConstraint : public IFunction
{
public:
	virtual double GetLimit () const = 0;
	virtual void SetLimit ( double rValue ) = 0;
//	virtual double Evaluate (const arma::vec& vCoordinates) const = 0;
	virtual bool IsSatisfied (const arma::vec& vSolution) const = 0;

	virtual EConstraintType GetType() const = 0;
	virtual ~IConstraint () {};
};

//////////////////////////////////////////////////

class FUNCTION ILinearConstraint : public IConstraint
{
public:
	virtual arma::vec GetCoeff () const = 0;
	virtual void SetCoeff (arma::vec& vCoef) = 0;
};

//////////////////////////////////////////////////