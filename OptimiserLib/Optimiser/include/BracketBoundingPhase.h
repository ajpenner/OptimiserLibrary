// BracketBoundingPhase.h
//
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "optimiser.h"
#include "optimisable.h"
#include <limits>
#include <armadillo>

//////////////////////////////////////////////////

namespace eBoundingPhaseParameters
{
	enum type
	{
		eDelta = 0
	};
}

class OPTIMISER CBracketBoundingPhase : public BBracket
{
public:
	CBracketBoundingPhase () {}
	~CBracketBoundingPhase () {}

	void SetLeftBracket (const arma::vec& v) override { m_LeftBracket = v; }
	void SetRightBracket (const arma::vec& v) override { m_RightBracket = v; }

	void SetParameters (const std::vector<variant>& parameters) override;

	void SetOptimisableFunction (std::shared_ptr<const IFunction> ptrOptimisableFunction) override;

	void Optimise (double /*target*/) override;

	arma::vec GetVector() const override { return arma::vec(); }
	arma::vec GetRange () const override;

private:
	arma::vec m_LeftBracket;
	arma::vec m_RightBracket;

	double m_delta;
};

//////////////////////////////////////////////////