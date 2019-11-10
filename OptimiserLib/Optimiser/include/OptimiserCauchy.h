// OptimiserCauchy.h
//
// 2016
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "optimiser.h"
#include "optimisable.h"
#include "LinearFunction.h"
#include "MinimiserGolden.h"

#include <limits>
#include <armadillo>

//////////////////////////////////////////////////

namespace eCauchyParameters
{
	enum type
	{
		eAlpha = 0
	};
}

class OPTIMISER COptimiserCauchy : public BOptimiser
{
public:
	COptimiserCauchy () : m_alpha(1.0) {}
	~COptimiserCauchy () {}

	void SetLeftBracket (const arma::vec& vLeft) override;
	void SetRightBracket (const arma::vec& vRight) override;

	void SetParameters (const std::vector<variant>& parameters) override;

	void Optimise (double target = 0) override;

	double GetInverseOOF () const override;

	void InverseOptimise (double target) override;

	double GetMinimalScaleValue ( std::shared_ptr<CLinearFunction> vf, CMinimiserGolden& minimiser, const arma::vec& gradient);

private:

	void CheckAndSwap (double& rLeft, double& rRight);

	double m_InverseOOFValue;

	arma::vec m_LeftBracket;
	arma::vec m_RightBracket;
	double m_alpha;
};

//////////////////////////////////////////////////