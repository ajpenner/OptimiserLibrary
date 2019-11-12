// OptimiserBox.cpp
//
// 2016
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif

//////////////////////////////////////////////////

#include "OptimiserBox.h"
#include "optimisable.h"
#include <cassert>

//////////////////////////////////////////////////

void COptimiserBox::SetLeftBracket (const arma::vec& vLeft)
{
	m_LeftBracket = vLeft; 
}

//////////////////////////////////////////////////

void COptimiserBox::SetRightBracket (const arma::vec& vRight) 
{ 
	m_RightBracket = vRight; 
}

//////////////////////////////////////////////////

void COptimiserBox::SetOptimisableFunction (std::shared_ptr<const IFunction> ptrOptimisable)
{
	m_ptrOptimisableFunction = ptrOptimisable;
}

//////////////////////////////////////////////////

double COptimiserBox::GetInverseOOF() const
{
	return m_InverseOOFValue;
}

//////////////////////////////////////////////////

void COptimiserBox::InverseOptimise(double target)
{
	Optimise(-target);
}

//////////////////////////////////////////////////

void COptimiserBox::Optimise (double target)
{
	assert (m_ptrOptimisableFunction);

	auto coordinates = std::vector<arma::vec> ();
	const arma::uword xIndex (0);
	const arma::uword yIndex (1);
	coordinates.emplace_back (arma::vec({ m_LeftBracket[xIndex], m_LeftBracket[yIndex] }));
	coordinates.emplace_back (arma::vec({ m_LeftBracket[xIndex], m_RightBracket[yIndex] }));
	coordinates.emplace_back (arma::vec({ m_RightBracket[xIndex], m_LeftBracket[yIndex] }));
	coordinates.emplace_back (arma::vec({ m_RightBracket[xIndex], m_RightBracket[yIndex] })); // provide reference to m_MinVector?

	// calculate the centre
	auto xC = (m_LeftBracket[xIndex] + m_RightBracket[xIndex])*0.5;
	auto yC = (m_LeftBracket[yIndex] + m_RightBracket[yIndex])*0.5;
	coordinates.emplace_back (arma::vec({ xC, yC }));

	auto delta = arma::vec({ m_LeftBracket[xIndex] - xC, m_LeftBracket[yIndex] - yC });

	m_OOFValue = m_ptrOptimisableFunction->Evaluate (coordinates[4]);

	// evaluate the box coordinates
	uint8_t index = 4; // assume the centre coordinate is correct, this will be fixed over time
	do
	{
		for (auto it = ++coordinates.crbegin(); it != coordinates.crend (); ++it)
		{
			auto f = m_ptrOptimisableFunction->Evaluate (*it);
			if (f < m_OOFValue)
			{
				m_OOFValue = f;
				index = static_cast<uint8_t>(std::distance (it, coordinates.crend ()) - 1);
			}
		}
		// now we have a minimal value and the corresponding index
		if (index == 4)
		{
			// nothing changed, we need to reduce the bounding box
			delta = delta / 2;
		}

		// define new coordinates using original delta
		coordinates[4] = coordinates[index];
		coordinates[0] = arma::vec ({ coordinates[4][0] - delta[0], coordinates[4][1] - delta[1] });
		coordinates[1] = arma::vec ({ coordinates[4][0] - delta[0], coordinates[4][1] + delta[1] });
		coordinates[2] = arma::vec ({ coordinates[4][0] + delta[0], coordinates[4][1] + delta[1] });
		coordinates[3] = arma::vec ({ coordinates[4][0] + delta[0], coordinates[4][1] - delta[1] });
		index = 4;
	} while (std::abs(m_OOFValue-target) > m_rTolerance);

	m_MinVector = coordinates[index];

	return;
}

//////////////////////////////////////////////////
