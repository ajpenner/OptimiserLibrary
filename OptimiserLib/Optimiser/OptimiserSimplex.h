// COptimiserSimplex.h
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

namespace eSimplexParameters
{
	enum type
	{
		eBeta = 0,
		eGamma,
		eSimplexStart
	};
}

class OPTIMISER COptimiserSimplex : public BOptimiser
{
public:
	COptimiserSimplex ();
	~COptimiserSimplex () {}

	void SetLeftBracket (const arma::vec& vLeft);
	void SetRightBracket (const arma::vec& vRight);

	void SetParameters (const std::vector<variant>& parameters) override;

	void Optimise (double target = 0) override;

	double GetInverseOOF () const override;

	void InverseOptimise (double target) override;

private:

	void CreateSimplex();
	
	double m_gamma;
	double m_beta;
	double m_delta;

	std::vector<arma::vec> m_simplex;

	double m_InverseOOFValue;

	arma::vec m_LeftBracket;
	arma::vec m_RightBracket;
};

//////////////////////////////////////////////////