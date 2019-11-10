// MinimiserExhaustive.cpp
//
// 2016
//////////////////////////////////////////////////

#include "stdafx.h"

//////////////////////////////////////////////////

#include "MinimiserExhaustive.h"
#include "optimisable.h"
#include <cassert>
#include <iostream>

//////////////////////////////////////////////////

void CMinimiserExhaustive::CheckAndSwap (double& rLeft, double& rRight)
{
	if (rRight < rLeft)
	{
		std::swap (rLeft, rRight);
		std::swap (m_LeftBracket, m_RightBracket);
	}
}

//////////////////////////////////////////////////

void CMinimiserExhaustive::Optimise(double /*target*/)
{
	assert (m_ptrOptimisableFunction);
	// evaluate end points, compare them, then use bisection to find the minimum
	auto delta = (m_RightBracket - m_LeftBracket) / static_cast<double>(m_iIntermediatePoints);

	auto rInc1 = (m_LeftBracket + delta).eval();
	auto rInc2 = (rInc1 + delta).eval();
	/*
	while (rInc2 != m_RightBracket) // need to sort this out to make this work. (I stopped because this method will never be used)
	{
		auto rValue1 = m_ptrOptimisableFunction->Evaluate (m_LeftBracket);
		auto rValue2 = m_ptrOptimisableFunction->Evaluate (rInc1);
		auto rValue3 = m_ptrOptimisableFunction->Evaluate (rInc2);

		if (rValue1 >= rValue2 && rValue3 >= rValue2)
		{
			m_DirectionVector = rInc1;
			return; // we are done
		}
		else
		{
			m_LeftBracket = rInc1.eval ();
			rInc1 = rInc2;
			rInc2 = rInc1 + delta;
		}
	}
	*/

	return;
}

//////////////////////////////////////////////////