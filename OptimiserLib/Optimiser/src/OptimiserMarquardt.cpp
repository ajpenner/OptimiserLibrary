// OptimiserMarquardt.cpp
//
// 2017
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif

//////////////////////////////////////////////////

#include "OptimiserMarquardt.h"
#include "optimisable.h"

#include "BracketBoundingPhase.h"

#include <cassert>

//////////////////////////////////////////////////

void COptimiserMarquardt::SetParameters (const std::vector<variant>& parameters)
{
	assert (parameters.size () >= 1);
	m_penalty = boost::get<double> (parameters[eMarquardtParameters::ePenalty]);
}

//////////////////////////////////////////////////

double COptimiserMarquardt::GetMinimalScaleValue(std::shared_ptr<CLinearFunction> ptrVF, CMinimiserGolden& minimiser, const arma::vec& direction)
{
	arma::vec guess = arma::vec({1});
	ptrVF->SetVector (direction);
	ptrVF->SetOffset (m_MinVector);

	minimiser.SetOptimisableFunction (ptrVF);
	
	CBracketBoundingPhase bp;
	std::vector<variant> delta = { 0.5 };
	bp.SetParameters (delta);
	bp.SetOptimisableFunction (ptrVF);
	bp.SetInitialGuess (guess);
	bp.Optimise (0);
	auto range = bp.GetRange ();
	
	minimiser.SetInitialGuess (guess);
	minimiser.SetLeftBracket (arma::vec ({ range.at(0) }));
	minimiser.SetRightBracket (arma::vec ({ range.at(1) }));

	minimiser.Optimise (0);

	return minimiser.GetVector ().at (0);
}

//////////////////////////////////////////////////

void COptimiserMarquardt::Optimise (double target)
{
	assert (m_ptrOptimisableFunction);

	auto vGrad = m_ptrOptimisableFunction->CalculateGradient(m_MinVector);
	if (arma::norm (vGrad) < m_rTolerance)
		return; // we are done

	auto ptrVF = std::make_shared<CLinearFunction>(m_ptrOptimisableFunction);
	CMinimiserGolden minimiser;
	minimiser.SetTolerance(m_rTolerance);

	arma::mat eye = arma::eye (m_ptrOptimisableFunction->GetFunctionDimension (), m_ptrOptimisableFunction->GetFunctionDimension ());
	auto mHess = m_ptrOptimisableFunction->CalculateHessian (m_MinVector);

	auto vNewVector = m_MinVector;

	uint32_t k = 0;
	auto rLambda = m_penalty;
	auto fx = m_ptrOptimisableFunction->Evaluate(m_MinVector);
	do
	{
		mHess = m_ptrOptimisableFunction->CalculateHessian (m_MinVector);
		vGrad = m_ptrOptimisableFunction->CalculateGradient (m_MinVector);

		arma::vec s = -(arma::inv_sympd(mHess + rLambda*eye)*vGrad).eval();
		rLambda = GetMinimalScaleValue (ptrVF, minimiser, s); // this provided a big improvement early in the iterations
		m_MinVector += rLambda*s;
		auto fxp1 = m_ptrOptimisableFunction->Evaluate(m_MinVector);
		if (fxp1 < fx)
		{
			m_OOFValue = fxp1;
			fx = fxp1;
			rLambda *= 0.5;
			++k;
		}
		else
		{
			rLambda *= 2;
		}
	} while (std::abs(m_OOFValue - target) > m_rTolerance && arma::norm(vGrad)>m_rTolerance);


	return;
}

//////////////////////////////////////////////////
