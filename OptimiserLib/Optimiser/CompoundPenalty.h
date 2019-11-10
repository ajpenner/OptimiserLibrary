// Compound.h
//
// A utility function that allows one to make compound functions
// This is useful in the penalty optimisation methods
//
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "optimisable.h"
#include "constraints.h"

#include <memory>
#include <cassert>

//////////////////////////////////////////////////

class FUNCTION CCompoundPenaltyFunction : public IFunction
{
public:
	CCompoundPenaltyFunction (std::shared_ptr<const IFunction> ptrOOF, std::vector<std::shared_ptr<const IConstraint>> vConstraints, std::vector<double> vCoeff);
	~CCompoundPenaltyFunction () {}

	double Evaluate (const arma::vec& vValues) const override;
	arma::vec CalculateGradient (const arma::vec& vec) const override;

	arma::mat CalculateHessian (const arma::vec& vec) const override;

	void SetCoeff (std::vector<double> vCoeff)
	{
		m_Coeff = vCoeff;
	}

	arma::uword GetFunctionDimension () const override;

private:
	std::shared_ptr<const IFunction> m_ptrOOF;
	std::vector<std::shared_ptr<const IConstraint>> m_Constraints;
	std::vector<double> m_Coeff;
};

//////////////////////////////////////////////////