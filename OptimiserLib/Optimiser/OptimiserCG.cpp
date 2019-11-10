// OptimiserCG.cpp
//
// 2017
//////////////////////////////////////////////////

#include "stdafx.h"

//////////////////////////////////////////////////

#include "OptimiserCG.h"
#include "optimisable.h"

#include "BracketBoundingPhase.h"

#include <cassert>
#include <iostream>

//////////////////////////////////////////////////

double COptimiserCG::GetMinimalScaleValue(std::shared_ptr<CLinearFunction> vf, CMinimiserGolden& minimiser, const arma::vec& direction)
{
	arma::vec guess = arma::vec({1});
	vf->SetVector (direction);
	vf->SetOffset (m_MinVector);

	minimiser.SetOptimisableFunction (vf);
	
	auto range = arma::vec ({0,1});
	
	minimiser.SetInitialGuess (guess);
	minimiser.SetLeftBracket (arma::vec ({ range.at(0) }));
	minimiser.SetRightBracket (arma::vec ({ range.at(1) }));

	minimiser.Optimise (0);

	return minimiser.GetVector ().at (0);
}

//////////////////////////////////////////////////

void COptimiserCG::Optimise (double target)
{
	assert (m_ptrOptimisableFunction);

	auto vGrad = m_ptrOptimisableFunction->CalculateGradient(m_MinVector);
	if (arma::norm (vGrad) < m_rTolerance)
		return; // we are done

	auto vf = std::make_shared<CLinearFunction>(m_ptrOptimisableFunction);
	CMinimiserGolden minimiser;
	minimiser.SetTolerance(m_rTolerance);

	uint32_t k = 0;
	auto rLambda = 1.0;
	double differential = realZero;
	auto fx = m_ptrOptimisableFunction->Evaluate(m_MinVector);
	vGrad = m_ptrOptimisableFunction->CalculateGradient (m_MinVector);
	arma::vec s = -vGrad;
	arma::vec vGradOriginal;
	do
	{
		rLambda = GetMinimalScaleValue (vf, minimiser, s);
		m_MinVector += rLambda*s;

		arma::vec originalPosition = m_MinVector; // need this to calculate the termination condition
		arma::vec normOriginal = arma::normalise(s);

		auto vGradp1 = m_ptrOptimisableFunction->CalculateGradient (m_MinVector);
		double vGradNorm = std::pow ( arma::norm (vGrad), 2);
		double vGradNormp1 = std::pow (arma::norm (vGradp1), 2);

		s = -vGradp1 + vGradNormp1 / vGradNorm *s;
		rLambda = GetMinimalScaleValue (vf, minimiser, s);
		arma::vec norm1 = arma::normalise(s);

		double linIndepMetric = arma::dot (normOriginal, norm1);
		if (linIndepMetric > 0.95) // less than ~18 degrees
		{
			// trigger restart
			s = -vGrad;
			++k;
			continue;
		}

		m_MinVector += rLambda*s;

		m_OOFValue = m_ptrOptimisableFunction->Evaluate(m_MinVector);

		differential = arma::norm (m_MinVector - originalPosition) / arma::norm (originalPosition);
		vGrad = m_ptrOptimisableFunction->CalculateGradient (m_MinVector);

		vGradNorm = std::pow (arma::norm (vGrad), 2);
		s = -vGrad + vGradNorm / vGradNormp1*s;
		++k;

	} while ( (differential > m_rTolerance && arma::norm(vGrad)>m_rTolerance) && (m_OOFValue-target) > m_rTolerance );

	return;
}

//////////////////////////////////////////////////