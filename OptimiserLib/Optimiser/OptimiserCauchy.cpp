// OptimiserCauchy.cpp
//
// 2017
//////////////////////////////////////////////////

#include "stdafx.h"

//////////////////////////////////////////////////

#include "OptimiserCauchy.h"
#include "optimisable.h"
#include "BracketBoundingPhase.h"

#include <cassert>
#include <iostream>

//////////////////////////////////////////////////

void COptimiserCauchy::SetLeftBracket (const arma::vec& vLeft)
{
	m_LeftBracket = vLeft; 
}

//////////////////////////////////////////////////

void COptimiserCauchy::SetRightBracket (const arma::vec& vRight)
{ 
	m_RightBracket = vRight; 
}

//////////////////////////////////////////////////

double COptimiserCauchy::GetInverseOOF() const
{
	return m_InverseOOFValue;
}

//////////////////////////////////////////////////

void COptimiserCauchy::SetParameters (const std::vector<variant>& parameters)
{
	assert (parameters.size () >= 1);
	m_alpha = boost::get<double>(parameters[eCauchyParameters::eAlpha]);
}

//////////////////////////////////////////////////

void COptimiserCauchy::CheckAndSwap(double& rLeft, double& rRight)
{
	if (std::abs (rLeft) < std::abs (rRight))
	{
		auto vTmp = m_LeftBracket;
		auto rTmp = rLeft;

		m_LeftBracket = m_RightBracket;
		rLeft = rRight;

		m_RightBracket = vTmp;
		rRight = rTmp;
	}
}

//////////////////////////////////////////////////

void COptimiserCauchy::InverseOptimise(double target)
{
	Optimise(-target);
}

//////////////////////////////////////////////////

double COptimiserCauchy::GetMinimalScaleValue ( std::shared_ptr<CLinearFunction> vf, CMinimiserGolden& minimiser, const arma::vec& gradient )
{
	arma::vec guess = -gradient;
	vf->SetVector(guess); // CLinearFunction needs to be renamed, perhaps CLinearFunction?
	vf->SetOffset(m_MinVector);

	minimiser.SetOptimisableFunction (vf);
/*
    CBracketBoundingPhase bp;
	bp.SetDelta (0.5);
	bp.SetOptimisableFunction (&vf);
	bp.SetInitialGuess (guess);
	bp.Optimise (0);
	auto thing = bp.GetRange ();
*/
	minimiser.SetInitialGuess(guess);
	minimiser.SetLeftBracket (arma::vec ({ 0 }));
	minimiser.SetRightBracket (arma::vec ({ 1 }));

	minimiser.Optimise(0);

	return minimiser.GetVector().at(0);
}

//////////////////////////////////////////////////

void COptimiserCauchy::Optimise( double target ) // incomplete
{
	assert (m_ptrOptimisableFunction);

	// Set minimiser
	CMinimiserGolden minimiser;
	auto ptrVF = std::make_shared<CLinearFunction>(m_ptrOptimisableFunction);
	minimiser.SetTolerance (1e-6);

	auto original = m_MinVector;
	uint32_t k = 0;

	do
	{
		auto gradient = m_ptrOptimisableFunction->CalculateGradient (m_MinVector);

		if (arma::norm (gradient) < m_rTolerance)
		{
			m_OOFValue = m_ptrOptimisableFunction->Evaluate (m_MinVector);
			return;
		}
		else if (k > m_MaxIteration) // give up 
		{
			m_OOFValue = m_ptrOptimisableFunction->Evaluate (m_MinVector);
			return;
		}
		else
		{
			original = m_MinVector;
			auto scale = GetMinimalScaleValue(ptrVF, minimiser, gradient);
			m_MinVector = m_MinVector - scale*gradient;
			++k;
		}
	} while (std::abs(arma::norm(m_MinVector - original)-target) > m_rTolerance);

	m_OOFValue = m_ptrOptimisableFunction->Evaluate (m_MinVector);

	return;
}

//////////////////////////////////////////////////