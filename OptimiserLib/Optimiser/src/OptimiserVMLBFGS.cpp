// OptimiserVMLBFGS.cpp
//
// 2016
//////////////////////////////////////////////////

#include "stdafx.h"

//////////////////////////////////////////////////

#include "OptimiserVMLBFGS.h"
#include "optimisable.h"

#include "LinearFunction.h"
#include "MinimiserGolden.h"
#include "BracketBoundingPhase.h"

#include "STLExtension.h"

#include <cassert>
#include <iostream>

//////////////////////////////////////////////////

void COptimiserVMLBFGS::SetParameters (const std::vector<variant>& parameters)
{
	m_Count = boost::get<int> (parameters[eLBFGSParameters::eCount]);
}

//////////////////////////////////////////////////

double COptimiserVMLBFGS::GetMinimalScaleValue (std::shared_ptr<CLinearFunction> ptrVF, CMinimiserGolden& minimiser, const arma::vec& direction)
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

arma::mat COptimiserVMLBFGS::UpdateHessian ()
{
	auto gamma = (m_deltaGrad.t()*m_deltaGrad).eval().at(0);
	assert (gamma != 0.0);
	m_mHess = m_deltaX.t ()*m_deltaGrad / gamma; // initial guess

	auto zBegin = boost::make_zip_iterator(boost::make_tuple(m_vDeltaGrad.crbegin (), m_vDeltaX.crbegin ()));
	auto zEnd = boost::make_zip_iterator(boost::make_tuple (m_vDeltaGrad.crend (), m_vDeltaX.crend()));

	arma::mat mV = m_mHess;
	for (auto it = zBegin; it != zEnd; ++it)
	{
		auto rDenom = (arma::trans (m_deltaGrad)*m_deltaX).eval ()(0);
		assert (rDenom != 0.0);
		auto rho = 1.0 / rDenom;

		auto& dGrad = it->get<eLBFGSValues::eDeltaGrad>();
		auto& dX = it->get<eLBFGSValues::eDeltaX>();
		mV *= (m_mEye - rho*dGrad*dX.t ());
	}

	for (auto it = zBegin; it != zEnd; ++it)
	{
		auto rDenom = (arma::trans (m_deltaGrad)*m_deltaX).eval ()(0);
		assert (rDenom != 0.0);
		auto rho = 1.0 / rDenom;

		auto& dGrad = it->get<eLBFGSValues::eDeltaGrad> ();
		auto& dX = it->get<eLBFGSValues::eDeltaX> ();
		mV *= (m_mEye - rho*dGrad*dX.t ());
	}
	return arma::mat (); // temporary
}

//////////////////////////////////////////////////

void COptimiserVMLBFGS::Optimise (double target)
{
	assert (m_ptrOptimisableFunction);

	const auto iDim = m_ptrOptimisableFunction->GetFunctionDimension ();
	m_mEye = arma::eye (iDim, iDim);

	m_vDeltaX.reserve(m_Count);
	m_vDeltaGrad.reserve(m_Count);

	auto vf = std::make_shared<CLinearFunction> (m_ptrOptimisableFunction);
	CMinimiserGolden minimiser;
	minimiser.SetTolerance (m_rTolerance);
	auto vGrad = m_ptrOptimisableFunction->CalculateGradient (m_MinVector);
	do
	{
		arma::vec s = -arma::normalise (m_mHess*vGrad);
		auto rLambda = GetMinimalScaleValue (vf, minimiser, s);

		if (m_vDeltaX.size () == m_Count)
		{
			pop_front(m_vDeltaX);
			pop_front(m_vDeltaGrad);
		}

		m_vDeltaX.push_back( rLambda*s );
		arma::vec newMin = m_MinVector + m_deltaX;
		arma::vec vGradp1 = m_ptrOptimisableFunction->CalculateGradient (newMin);
		m_vDeltaGrad.push_back(vGradp1 - vGrad);

		m_mHess = UpdateHessian();

		vGrad = vGradp1;
		m_MinVector = newMin;
		m_OOFValue = m_ptrOptimisableFunction->Evaluate (m_MinVector);
	} while (arma::norm (vGrad) > m_rTolerance && m_OOFValue > m_rTolerance);

	return;
}

//////////////////////////////////////////////////