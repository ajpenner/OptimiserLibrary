// OptimiserBox.h
//
// 2016
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "optimiser.h"
#include "optimisable.h"
#include <limits>
#include <armadillo>

//////////////////////////////////////////////////

class OPTIMISER COptimiserBox : public BOptimiser
{
public:
	COptimiserBox () {}
	~COptimiserBox () {}

	void SetLeftBracket (const arma::vec& vLeft) override;
	void SetRightBracket (const arma::vec& vRight) override;

	void SetParameters (const std::vector<variant>& parameters) override {};

	void SetOptimisableFunction (std::shared_ptr<const IFunction> ptrOptimisableFunction) override;

	void Optimise (double target = 0) override;

	double GetInverseOOF () const override;

	void InverseOptimise (double target) override;

private:

	double m_InverseOOFValue;

	arma::vec m_LeftBracket;
	arma::vec m_RightBracket;
};

//////////////////////////////////////////////////