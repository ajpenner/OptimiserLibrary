// MinimiserGolden.cpp
//
// 2016
//////////////////////////////////////////////////

#include "stdafx.h"

//////////////////////////////////////////////////

#include "MinimiserGolden.h"
#include "optimisable.h"
#include <cassert>
#include <iostream>

//////////////////////////////////////////////////

void CMinimiserGolden::CalculateUpperBracket (arma::vec& vUpperBracket)
{
	vUpperBracket = (m_RightBracket - m_rGoldenRatio *(m_RightBracket - m_LeftBracket)).eval ();
}

//////////////////////////////////////////////////

void CMinimiserGolden::CalculateLowerBracket (arma::vec& vLowerBracket)
{
	vLowerBracket = (m_LeftBracket + m_rGoldenRatio*(m_RightBracket - m_LeftBracket)).eval ();
}

//////////////////////////////////////////////////

void CMinimiserGolden::Optimise (double /*target*/)
{
	assert (m_ptrOptimisableFunction);

	Transform transformed(m_ptrOptimisableFunction, m_LeftBracket.at(0), m_RightBracket.at(0));
	// the transformed function is the one we are optimising below
	m_LeftBracket = arma::vec({0});
	m_RightBracket = arma::vec({1});
	auto rangeWidth = 1.0;

	// need inverse transform to get the final answer

	// evaluate end points, compare them, then use bisection to find the minimum
	auto w1 = m_LeftBracket.at(0) + m_rGoldenRatio * rangeWidth;
	auto f1 = transformed.Evaluate (arma::vec ({ w1 }));
	auto w2 = m_RightBracket.at(0) - m_rGoldenRatio * rangeWidth;
	auto f2 = transformed.Evaluate(arma::vec({ w2 }));

	do
	{
		if ( f1 >= f2 )
		{
			m_RightBracket = w1;
			rangeWidth = m_RightBracket.at(0) - m_LeftBracket.at(0);
			w1 = w2;
			f1 = f2;

			w2 = m_RightBracket.at(0) - m_rGoldenRatio*rangeWidth;
			f2 = transformed.Evaluate (arma::vec ({ w2 }));
		}
		else
		{
			m_LeftBracket = w2;

			rangeWidth = m_RightBracket.at(0) - m_LeftBracket.at(0);
			w2 = w1;	
			f2 = f1;

			w1 = m_LeftBracket.at(0) + m_rGoldenRatio*rangeWidth;
			f1 = transformed.Evaluate (arma::vec ({w1}));
		}

	} while ( rangeWidth > m_rTolerance);

	// may as well average the w's and get the original coordinate back
	auto aveW = 0.5*(w1 + w2);
	m_MinVector = transformed.GetOriginalCoordinate (aveW);
	m_OOFValue = m_ptrOptimisableFunction->Evaluate (m_MinVector);

	return;
}

//////////////////////////////////////////////////