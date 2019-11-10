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

namespace eMOMCompoundParameters
{
	enum type
	{
		eConstraint,
		eMultiplier
	};
}

class FUNCTION CCompoundMOMFunction : public IFunction
{
public:
	CCompoundMOMFunction (std::shared_ptr<const IFunction> ptrOOF, std::vector<std::shared_ptr<const IConstraint>> vConstraints, std::vector<double> vCoeff);
	~CCompoundMOMFunction () {}

	double Evaluate (const arma::vec& vValues) const override;
	arma::vec CalculateGradient (const arma::vec& vec) const override;
	arma::mat CalculateHessian (const arma::vec& vec) const override;
	arma::uword GetFunctionDimension () const override;

	void SetMultipliers (const arma::vec& vec);

	std::vector<double> GetConstraintViolations (const arma::vec& vec) const;

	void SetPenaltyFactor (double penalty);

private:
	std::shared_ptr<const IFunction> m_ptrOOF;
	std::vector<std::shared_ptr<const IConstraint>> m_Constraints;
	std::vector<double> m_Multiplier;
	double m_rPenalty;
};

//////////////////////////////////////////////////