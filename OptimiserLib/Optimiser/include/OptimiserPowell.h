// OptimiserPowell.h
//
// 2016
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "optimiser.h"
#include "optimisable.h"
#include <limits>

class CLinearFunction;
class CMinimiserGolden;
class CBracketBoundingPhase;

//////////////////////////////////////////////////

namespace ePowellParameters
{
	enum type
	{
		eAlpha = 0
	};
}

class OPTIMISER COptimiserPowell : public BOptimiser
{
public:
	COptimiserPowell() : m_alpha(1.0) {}
	~COptimiserPowell() {}

	void SetLeftBracket (const arma::vec& /*vec*/) override { return; }
	void SetRightBracket (const arma::vec& /*vec*/) override { return; }

	void SetParameters (const std::vector<variant>& parameters) override;

	void Optimise(double target = 0) override;

	double GetMinimalScaleValue (std::shared_ptr<CLinearFunction> ptrVF, CMinimiserGolden& minimiser, CBracketBoundingPhase& range, const arma::vec& direction);

private:

	double m_alpha;
};

//////////////////////////////////////////////////