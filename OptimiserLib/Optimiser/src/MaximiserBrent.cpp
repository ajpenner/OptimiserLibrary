// MinimiserBrent.cpp
//
// 2017
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif

//////////////////////////////////////////////////

#include "MaximiserBrent.h"
#include "optimisable.h"
#include <cassert>
#include <iostream>

//////////////////////////////////////////////////

void CMaximiserBrent::_CheckAndSwap (double& rLeft, double& rRight)
{
	if (rRight < rLeft)
	{
		auto vTmp = m_LeftBracket;
		auto rTmp = rLeft;

		m_LeftBracket = m_RightBracket;
		rLeft = rRight;

		m_RightBracket = vTmp;
		rRight = rTmp;
	}
}

//////////////////////////////////////////////////

void CMaximiserBrent::_PerformSearch(double& rLeft, double& rRight)
{
	auto vMiddleBracket = m_LeftBracket;
	auto rMid = -m_ptrOptimisableFunction->Evaluate (vMiddleBracket);
	auto bFlag = true;

	auto vBrent = arma::vec (m_LeftBracket.size ());
	vBrent.fill (std::numeric_limits<double>::max());

	do
	{
		auto vValueBracket = arma::vec (m_LeftBracket.size ());
		vValueBracket.fill (std::numeric_limits<double>::max());
		if (rLeft != rMid && rRight != rMid)
		{
			vValueBracket = (m_LeftBracket* rRight*rMid / (rLeft - rRight) / (rLeft - rMid)).eval ();
			vValueBracket += (m_RightBracket*rLeft*rMid / (rLeft - rRight) / (rRight - rMid)).eval ();
			vValueBracket += (vMiddleBracket*rLeft*rRight / (rMid - rLeft) / (rMid - rRight)).eval ();
		}
		else
		{
			vValueBracket = m_RightBracket - rRight*(m_RightBracket - m_LeftBracket) / (rRight - rLeft); // secant
		}

		if (
			!((arma::norm (3.0*m_LeftBracket + m_RightBracket / 4.0) < arma::norm (vValueBracket)) && (arma::norm (vValueBracket) < arma::norm (m_RightBracket))) ||
			(bFlag && arma::norm (vValueBracket - m_RightBracket) >= arma::norm (m_RightBracket - vMiddleBracket) / 2.0) ||
			(!bFlag && arma::norm (vValueBracket - m_RightBracket) >= arma::norm (vMiddleBracket - vBrent) / 2.0) ||
			(bFlag && (arma::norm (m_RightBracket - vMiddleBracket) < m_rTolerance)) ||
			(!bFlag && arma::norm (vMiddleBracket - vBrent) < m_rTolerance)
			)
		{
			vValueBracket = (m_LeftBracket + m_RightBracket) / 2.0; // bisection
			bFlag = true;
		}
		else
		{
			bFlag = false;
		}

		auto rValue = -m_ptrOptimisableFunction->Evaluate (vValueBracket);
		vBrent = vMiddleBracket;
		vMiddleBracket = m_RightBracket;

		if (rLeft < rValue)
		{
			m_RightBracket = vValueBracket;
			rRight = -m_ptrOptimisableFunction->Evaluate (m_RightBracket);
		}
		else
		{
			m_LeftBracket = vValueBracket;
			rLeft = -m_ptrOptimisableFunction->Evaluate (m_LeftBracket);
		}

		_CheckAndSwap (rLeft, rRight);
		vMiddleBracket = m_LeftBracket;
		rMid = -m_ptrOptimisableFunction->Evaluate (vMiddleBracket);
	} while (std::abs (rLeft - rRight) > m_rTolerance && arma::norm (m_LeftBracket - m_RightBracket) > m_rTolerance);
}

//////////////////////////////////////////////////

void CMaximiserBrent::Optimise (double /*target*/)
{
	assert (m_ptrOptimisableFunction);
	// evaluate end points, compare them, then use bisection to find the minimum
	auto vOriginalRight = m_RightBracket;
	auto rLeft = -m_ptrOptimisableFunction->Evaluate(m_LeftBracket);
	auto rRight = -m_ptrOptimisableFunction->Evaluate(m_RightBracket);
	auto rOriginalRight = rRight;

	auto bIsSplit = false;

	auto vSplit = m_RightBracket;

	if (rRight == rLeft)
	{
		int iIsSame = 0;
		for (arma::uword i = 0; i < m_RightBracket.size (); ++i)
		{
			if (m_RightBracket[i] == m_LeftBracket[i])
			{
				++iIsSame;
			}
		}
		if (iIsSame == m_LeftBracket.size ())
		{
			std::cout << "Use a valid bracket range" << std::endl;
			return;
		}
		else // below we arbitrarily throw away half of the search space. Will need to evaluate both halves to ensure a user-expected minimum
		{
			m_RightBracket = 0.5*(m_LeftBracket + m_RightBracket).eval ();
			rRight = -m_ptrOptimisableFunction->Evaluate (m_RightBracket);
			bIsSplit = true;
		}
	}

	_CheckAndSwap (rLeft, rRight);

	_PerformSearch (rLeft, rRight);
	if (rLeft < rRight)
	{
		m_MinVector = m_LeftBracket;
		m_OOFValue = rLeft;
	}
	else
	{
		m_MinVector = m_RightBracket;
		m_OOFValue = rRight;
	}

	// did the first half, now, if necessary do second half
	if (bIsSplit == true)
	{
		rLeft = m_OOFValue;
		m_LeftBracket = m_MinVector;
		m_RightBracket = vOriginalRight;
		rRight = rOriginalRight;
		_CheckAndSwap(rLeft, rRight);
		_PerformSearch (rLeft, rRight);

		// these are reltive to the first half oof's
		if (rLeft < rRight)
		{
			m_MinVector = m_LeftBracket;
			m_OOFValue = rLeft;
		}
		else
		{
			m_MinVector = m_RightBracket;
			m_OOFValue = rRight;
		}
	}

	return;
}

//////////////////////////////////////////////////
