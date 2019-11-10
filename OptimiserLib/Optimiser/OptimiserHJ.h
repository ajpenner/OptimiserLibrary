// COptimiserHJ.h
// Hooke - Jeeves
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "optimiser.h"
#include "optimisable.h"
#include <limits>
#include <armadillo>

//////////////////////////////////////////////////

namespace eHJParameters
{
	enum type
	{
		eIncrement = 0,
		eReduction
	};
}

class OPTIMISER COptimiserHJ : public BOptimiser
{
public:
	COptimiserHJ ();
	~COptimiserHJ () {}

	void SetLeftBracket (const arma::vec& vLeft) override;
	void SetRightBracket (const arma::vec& vRight) override;

	void SetParameters (const std::vector<variant>& parameters) override;

	void Optimise (double target = 0) override;

	double GetInverseOOF () const override;

	void InverseOptimise (double target) override;

	arma::vec Explore(arma::vec& start, bool& bSuccess);

private:

	void CreateSimplex();
	
	double m_ReductionFactor;
	double m_IncrementFactor;

	std::vector<arma::vec> m_simplex;

	double m_InverseOOFValue;

	arma::vec m_LeftBracket;
	arma::vec m_RightBracket;

	arma::mat m_I;
};

//////////////////////////////////////////////////