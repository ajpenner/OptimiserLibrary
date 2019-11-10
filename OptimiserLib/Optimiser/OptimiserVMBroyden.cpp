// OptimiserVMBroyden.cpp
//
// 2017
//////////////////////////////////////////////////

#include "stdafx.h"

//////////////////////////////////////////////////

#include "OptimiserVMBroyden.h"
#include "optimisable.h"

#include "LinearFunction.h"
#include "MinimiserGolden.h"
#include "BracketBoundingPhase.h"

#include <cassert>
#include <iostream>

//////////////////////////////////////////////////

arma::mat COptimiserVMBroyden::UpdateHessian()
{
	auto rDenom = (arma::trans (m_deltaGrad)*m_deltaX).eval ()(0);
	assert (rDenom != 0.0);
	auto rho = 1.0 / rDenom;

	arma::mat mP2 = rho*m_deltaGrad*m_deltaGrad.t ();

	double scalar = (m_deltaX.t () * m_mHess * m_deltaX).eval().at(0);
	assert (scalar != 0);

	arma::vec Bs = m_mHess*m_deltaX;
	arma::mat mP1 = Bs*Bs.t()/scalar;
	arma::vec v = m_deltaGrad*rho - Bs/scalar;

	arma::mat mP3 = m_phi*scalar*v*v.t();

	return m_mHess - mP1 + mP2 + mP3;
}

//////////////////////////////////////////////////

void COptimiserVMBroyden::SetParameters (const std::vector<variant>& parameters)
{
	m_phi = boost::get<double> (parameters[eBroydenParameters::ePhi]);
	assert (m_phi >= realZero && m_phi <= 1.0);
}

double COptimiserVMBroyden::GetMinimalScaleValue (std::shared_ptr<CLinearFunction> ptrVF, CMinimiserGolden& minimiser, const arma::vec& direction)
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

void COptimiserVMBroyden::Optimise (double target)
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
		m_mInvHess = arma::inv_sympd (m_mHess);

		vGrad = vGradp1;
		m_MinVector = newMin;
		m_OOFValue = m_ptrOptimisableFunction->Evaluate (m_MinVector);
	} while (arma::norm (vGrad) > m_rTolerance && m_OOFValue > m_rTolerance);

	return;
}

//////////////////////////////////////////////////