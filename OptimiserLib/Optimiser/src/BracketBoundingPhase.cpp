// BracketBoundingPhase.cpp
//
// 2017
//////////////////////////////////////////////////

#include "stdafx.h"

//////////////////////////////////////////////////

#include "BracketBoundingPhase.h"
#include "optimisable.h"
#include <cassert>
#include <iostream>

//////////////////////////////////////////////////

void CBracketBoundingPhase::SetParameters (const std::vector<variant>& parameters)
{
	m_delta = boost::get<double> (parameters[eBoundingPhaseParameters::eDelta]);
}

//////////////////////////////////////////////////

void CBracketBoundingPhase::SetOptimisableFunction (std::shared_ptr<const IFunction> ptrOptimisable)
{
	m_ptrOptimisableFunction = ptrOptimisable;
}

//////////////////////////////////////////////////

arma::vec CBracketBoundingPhase::GetRange() const
{
	if (m_LeftBracket.at(0) < m_RightBracket.at(0))
	{
		return arma::vec ({ m_LeftBracket.at(0), m_RightBracket.at(0) });
	}
	else
	{
		return arma::vec ({ m_RightBracket.at(0), m_LeftBracket.at(0) });

	}
}

//////////////////////////////////////////////////

void CBracketBoundingPhase::Optimise (double /*target*/)
{
	assert (m_ptrOptimisableFunction);

	uint32_t k = 0;

	auto fx = m_ptrOptimisableFunction->Evaluate (m_MinVector);
	m_RightBracket = m_MinVector + std::abs (m_delta);
	auto fxp1 = m_ptrOptimisableFunction->Evaluate (m_RightBracket);
	m_LeftBracket = m_MinVector - std::abs (m_delta);
	auto fxm1 = m_ptrOptimisableFunction->Evaluate (m_LeftBracket);

	if (fxp1 <= fx && fx <= fxm1)
	{
		m_delta = std::abs(m_delta);
	}
	else if(fxm1 <= fx && fx <= fxp1)
	{
		m_delta = -std::abs(m_delta);
	}
	else if( fxm1 > fx && fx < fxp1)
	{
		return;
	}

	do
	{
		m_RightBracket = m_MinVector + std::pow(2, k)*m_delta;
		fxp1 = m_ptrOptimisableFunction->Evaluate(m_RightBracket);

		if (fxp1 < fx)
		{
			++k;
			fx = fxp1;
			m_LeftBracket = m_MinVector;
			m_MinVector = m_RightBracket;
		}
		else
		{
			break;
		}
	} while (true); // forever

	return;
}

//////////////////////////////////////////////////