// OptimiserVMBFGS.cpp
//
// 2016
//////////////////////////////////////////////////

#include "stdafx.h"

//////////////////////////////////////////////////

#include "OptimiserVMSR1.h"
#include "optimisable.h"

#include "LinearFunction.h"
#include "BracketBoundingPhase.h"
#include "MinimiserGolden.h"

#include <cassert>
#include <iostream>

//////////////////////////////////////////////////

arma::mat COptimiserVMSR1::UpdateHessian ()
{
	// prepare new Hessian approx
	arma::vec v = (m_deltaX - m_mHess*m_deltaGrad).eval();
	double rDenom = (v.t()*m_deltaGrad).eval ()(0);
	assert (rDenom != 0.0);

	return (v*v.t()).eval() / rDenom;
}

//////////////////////////////////////////////////

double COptimiserVMSR1::GetMinimalScaleValue (std::shared_ptr<CLinearFunction> ptrVF, CMinimiserGolden& minimiser, const arma::vec& direction)
{
	arma::vec guess = arma::vec ({ 1 });
	ptrVF->SetVector (direction);
	ptrVF->SetOffset (m_MinVector);

	minimiser.SetOptimisableFunction (ptrVF);

	CBracketBoundingPhase bp;
	std::vector<variant> delta = { 0.5 };
	bp.SetParameters (delta);
	bp.SetOptimisableFunction (ptrVF);
	bp.SetInitialGuess (arma::vec{ 0 });
	bp.Optimise (0);
	auto range = bp.GetRange ();

	minimiser.SetInitialGuess (guess);
	minimiser.SetLeftBracket (arma::vec ({ range.at (0) }));
	minimiser.SetRightBracket (arma::vec ({ range.at (1) }));

	minimiser.Optimise (0);

	return minimiser.GetVector ().at (0);
}

//////////////////////////////////////////////////

void COptimiserVMSR1::Optimise (double target)
{
	assert (m_ptrOptimisableFunction);

	const auto iDim = m_ptrOptimisableFunction->GetFunctionDimension();
	m_mHess = arma::eye(iDim, iDim);

	auto ptrVF = std::make_shared<CLinearFunction>(m_ptrOptimisableFunction);
	CMinimiserGolden minimiser;
	minimiser.SetTolerance (m_rTolerance);

	auto vGrad = m_ptrOptimisableFunction->CalculateGradient (m_MinVector);
	do
	{
		arma::vec s = -arma::normalise(m_mHess*vGrad);
		auto rLambda = GetMinimalScaleValue (ptrVF, minimiser, s);

		m_deltaX = rLambda*s;
		arma::vec newMin = m_MinVector + m_deltaX;
		arma::vec vGradp1 = m_ptrOptimisableFunction->CalculateGradient (newMin);
		m_deltaGrad = vGradp1 - vGrad;

		m_mHess += UpdateHessian();

		vGrad = vGradp1;
		m_MinVector = newMin;
		m_OOFValue = m_ptrOptimisableFunction->Evaluate (m_MinVector);
	} while (arma::norm (vGrad) > m_rTolerance && m_OOFValue > m_rTolerance);

	return;
}

//////////////////////////////////////////////////