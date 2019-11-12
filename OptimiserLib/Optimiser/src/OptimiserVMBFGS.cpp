// OptimiserVMBFGS.cpp
//
// 2016
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif

//////////////////////////////////////////////////

#include "OptimiserVMBFGS.h"
#include "optimisable.h"

#include "LinearFunction.h"
#include "MinimiserGolden.h"
#include "BracketBoundingPhase.h"

#include <cassert>

//////////////////////////////////////////////////

double COptimiserVMBFGS::GetMinimalScaleValue (std::shared_ptr<CLinearFunction> ptrVF, CMinimiserGolden& minimiser, const arma::vec& direction)
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

arma::mat COptimiserVMBFGS::UpdateHessian()
{
	auto rDenom = (arma::trans (m_deltaGrad)*m_deltaX).eval ()(0);
	assert (rDenom != 0.0);
	auto rho = 1.0 / rDenom;

	arma::mat mP3 = rho*m_deltaX*m_deltaX.t ();

	arma::mat mP1 = m_mEye - rho*m_deltaX*m_deltaGrad.t ();
	arma::mat mP2 = m_mEye - rho*m_deltaGrad*m_deltaX.t ();

	return mP1*m_mHess*mP2 + mP3;
}

//////////////////////////////////////////////////

void COptimiserVMBFGS::Optimise (double target)
{
	assert (m_ptrOptimisableFunction);

	const auto iDim = m_ptrOptimisableFunction->GetFunctionDimension();
	m_mEye = arma::eye(iDim, iDim);
	m_mHess = m_mEye;

	auto vf = std::make_shared<CLinearFunction> (m_ptrOptimisableFunction);
	CMinimiserGolden minimiser;
	minimiser.SetTolerance (m_rTolerance);
	auto vGrad = m_ptrOptimisableFunction->CalculateGradient (m_MinVector);
	do
	{
		arma::vec s = -arma::normalise(m_mHess*vGrad);
		auto rLambda = GetMinimalScaleValue (vf, minimiser, s);

		m_deltaX = rLambda*s;
		arma::vec newMin = m_MinVector + m_deltaX;
		arma::vec vGradp1 = m_ptrOptimisableFunction->CalculateGradient (newMin);
		m_deltaGrad = vGradp1 - vGrad;

		m_mHess = UpdateHessian();

		vGrad = vGradp1;
		m_MinVector = newMin;
		m_OOFValue = m_ptrOptimisableFunction->Evaluate (m_MinVector);
	} while (arma::norm (vGrad) > m_rTolerance && m_OOFValue > m_rTolerance);

	return;
}

//////////////////////////////////////////////////
