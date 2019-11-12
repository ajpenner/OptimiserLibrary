// MinimiserGolden.cpp
//
// 2016
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif

//////////////////////////////////////////////////

#include "MinimiserBisection.h"
#include "optimisable.h"
#include <cassert>

//////////////////////////////////////////////////

void CMinimiserBisection::Optimise (double /*target*/)
{
	assert (m_ptrOptimisableFunction);
	// evaluate end points, compare them, then use bisection to find the minimum
	do
	{
		auto rLeft = m_ptrOptimisableFunction->Evaluate(m_LeftBracket);
		auto rRight = m_ptrOptimisableFunction->Evaluate (m_RightBracket);

		arma::vec vMiddleBracket = 0.5*(m_LeftBracket + m_RightBracket).eval();
		auto rMiddle = m_ptrOptimisableFunction->Evaluate (vMiddleBracket);

		if ( rMiddle < rRight  && rLeft < rMiddle )
		{
			m_RightBracket = vMiddleBracket;
			rRight = rMiddle;
		}
		else
		{
			m_LeftBracket = vMiddleBracket;
			rLeft = rMiddle;
		}

	} while (arma::norm (m_LeftBracket - m_RightBracket) > m_rTolerance);

	if (arma::norm (m_LeftBracket) < arma::norm (m_RightBracket))
	{
		m_MinVector = m_LeftBracket;
	}
	else
	{
		m_MinVector = m_RightBracket;
	}

	m_OOFValue = m_ptrOptimisableFunction->Evaluate (m_MinVector);

	return;
}

//////////////////////////////////////////////////
