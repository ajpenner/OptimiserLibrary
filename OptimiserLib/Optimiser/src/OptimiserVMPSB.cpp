// OptimiserVMPSB.cpp
//
// 2017
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif

//////////////////////////////////////////////////

#include "OptimiserVMPSB.h"
#include "optimisable.h"

#include "LinearFunction.h"
#include "MinimiserGolden.h"
#include "BracketBoundingPhase.h"

#include <cassert>

//////////////////////////////////////////////////

arma::mat COptimiserVMPSB::UpdateHessian()
{
	auto rDenom = (m_deltaX.t()*m_deltaX).eval ()(0);
	assert (rDenom != 0.0);
	auto rho = 1.0 / rDenom;

	arma::vec V = m_deltaGrad - m_mHess*m_deltaX;
	arma::mat mP1 = V*m_deltaX.t()*rho;
	arma::mat mP2 = m_deltaX*V.t()*rho;

	arma::mat mP3 = m_deltaX*(m_deltaX.t ()*V)*m_deltaX.t()*rho*rho;

	return m_mHess + mP1 + mP2 - mP3;
}

//////////////////////////////////////////////////

void COptimiserVMPSB::SetParameters (const std::vector<variant>& parameters)
{
	return;
}

double COptimiserVMPSB::GetMinimalScaleValue (std::shared_ptr<CLinearFunction> ptrVF, CMinimiserGolden& minimiser, const arma::vec& direction)
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

void COptimiserVMPSB::Optimise (double target)
{
	assert (m_ptrOptimisableFunction);

	const auto iDim = m_ptrOptimisableFunction->GetFunctionDimension();
	m_mEye = arma::eye(iDim, iDim);
	m_mHess = m_mEye;
	arma::mat m_mInvHess = m_mEye;

	auto ptrVF = std::make_shared<CLinearFunction> (m_ptrOptimisableFunction);
	CMinimiserGolden minimiser;
	minimiser.SetTolerance (m_rTolerance);
	auto vGrad = m_ptrOptimisableFunction->CalculateGradient (m_MinVector);

	do
	{
		arma::vec s = arma::normalise (-m_mInvHess * vGrad);
		auto rLambda = GetMinimalScaleValue (ptrVF, minimiser, s);

		m_deltaX = rLambda*s;
		arma::vec newMin = m_MinVector + m_deltaX;
		arma::vec vGradp1 = m_ptrOptimisableFunction->CalculateGradient (newMin);
		m_deltaGrad = (vGradp1 - vGrad).eval ();

		m_mHess = UpdateHessian();
		m_mInvHess = m_mHess.i();

		vGrad = vGradp1;
		m_MinVector = newMin;
		m_OOFValue = m_ptrOptimisableFunction->Evaluate (m_MinVector);
	} while (arma::norm (vGrad) > m_rTolerance && m_OOFValue > m_rTolerance);

	return;
}

//////////////////////////////////////////////////
