// OptimiserNewton.cpp
//
// 2017
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif

//////////////////////////////////////////////////

#include "OptimiserNewton.h"
#include "optimisable.h"
#include <cassert>

//////////////////////////////////////////////////

void COptimiserNewton::Optimise (double target)
{
	assert (m_ptrOptimisableFunction);

	auto mHess = m_ptrOptimisableFunction->CalculateHessian(m_MinVector);
	auto vGrad = m_ptrOptimisableFunction->CalculateGradient(m_MinVector);

	auto vNewVector = m_MinVector;

	auto rAlpha = 1.0e-1; // this is not correct

	do
	{
		m_MinVector -= (arma::inv_sympd (mHess)*vGrad).eval();
		mHess = m_ptrOptimisableFunction->CalculateHessian (m_MinVector);
		vGrad = m_ptrOptimisableFunction->CalculateGradient (m_MinVector);
	} while ((m_OOFValue = m_ptrOptimisableFunction->Evaluate (m_MinVector)) - target > m_rTolerance);


	return;
}

//////////////////////////////////////////////////
