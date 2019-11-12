// OptimiserBrent.cpp
//
// 2016
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif

//////////////////////////////////////////////////

#include "OptimiserBrent.h"
#include "optimisable.h"
#include <cassert>

//////////////////////////////////////////////////

void COptimiserBrent::SetLeftBracket (const arma::vec& vLeft)
{
	m_LeftBracket = vLeft; 
}

//////////////////////////////////////////////////

void COptimiserBrent::SetRightBracket (const arma::vec& vRight) 
{ 
	m_RightBracket = vRight; 
}

//////////////////////////////////////////////////

double COptimiserBrent::GetInverseOOF() const
{
	return m_InverseOOFValue;
}

//////////////////////////////////////////////////

void COptimiserBrent::CheckAndSwap(double& rLeft, double& rRight)
{
	if (std::abs (rLeft) < std::abs (rRight))
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

void COptimiserBrent::InverseOptimise(double target)
{
	Optimise(-target);
}

//////////////////////////////////////////////////

void COptimiserBrent::Optimise( double target )
{
	assert (m_ptrOptimisableFunction);
	// evaluate end points, compare them, then use bisection to find the minimum
	auto rLeft = m_ptrOptimisableFunction->Evaluate (m_LeftBracket);
	auto rRight = m_ptrOptimisableFunction->Evaluate (m_RightBracket);
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
		else
		{
			m_RightBracket = 0.5*(m_LeftBracket + m_RightBracket).eval ();
			rRight = m_ptrOptimisableFunction->Evaluate (m_RightBracket);
		}
	}
	if (rLeft*rRight > 0)
	{
		std::cout << "Bracket does not contain root" << std::endl;
		return;
	}
	CheckAndSwap(rLeft, rRight);

	auto vMiddleBracket = m_LeftBracket;
	auto rMid = m_ptrOptimisableFunction->Evaluate (vMiddleBracket);
	auto bFlag = true;
	
	auto vBrent = arma::vec (m_LeftBracket.size ());
	vBrent.fill(std::numeric_limits<double>::max());

	do
	{
		auto vValueBracket = arma::vec (m_LeftBracket.size ());
		vValueBracket.fill ( std::numeric_limits<double>::max() );
		if (rLeft != rMid && rRight != rMid)
		{
			vValueBracket = (m_LeftBracket* rRight*rMid / (rLeft - rRight) / (rLeft - rMid)).eval();
			vValueBracket += (m_RightBracket*rLeft*rMid / (rLeft - rRight) / (rRight - rMid)).eval();
			vValueBracket += (vMiddleBracket*rLeft*rRight / (rMid - rLeft) / (rMid - rRight)).eval();
		}
		else
		{
			vValueBracket = m_RightBracket - rRight*(m_RightBracket - m_LeftBracket) / (rRight - rLeft); // secant
		}

		if ( // I think brackets and function eval got switched
			!( (arma::norm(3.0*m_LeftBracket + m_RightBracket / 4.0) < arma::norm(vValueBracket)) && ( arma::norm(vValueBracket) < arma::norm(m_RightBracket)) )||
			(bFlag && arma::norm(vValueBracket - m_RightBracket) >= arma::norm(m_RightBracket - vMiddleBracket) / 2.0) ||
			(!bFlag && arma::norm(vValueBracket - m_RightBracket) >= arma::norm(vMiddleBracket - vBrent) / 2.0) ||
			(bFlag && (arma::norm(m_RightBracket - vMiddleBracket) < m_rTolerance)) ||
			(!bFlag && arma::norm(vMiddleBracket - vBrent) < m_rTolerance)
			)
		{
			vValueBracket = (m_LeftBracket + m_RightBracket) / 2.0; // bisection
			bFlag = true;
		}
		else
		{
			bFlag = false;
		}

		auto rValue = m_ptrOptimisableFunction->Evaluate(vValueBracket);
		vBrent = vMiddleBracket;
		vMiddleBracket = rRight;

		if ( target * rLeft*rValue < 0)
		{
			m_RightBracket = vValueBracket;
			rRight = m_ptrOptimisableFunction->Evaluate (m_RightBracket);
		}
		else
		{
			m_LeftBracket = vValueBracket;
			rLeft = m_ptrOptimisableFunction->Evaluate (m_LeftBracket);
		}

//		CheckAndSwap (rLeft, rRight);

		if (std::abs(rLeft) < m_rTolerance)
			break;

		if (std::abs(rRight) < m_rTolerance)
			break;

	} while (std::abs(arma::norm (m_LeftBracket - m_RightBracket)-target) > m_rTolerance);

	if (std::abs (rLeft) < m_rTolerance && std::abs (rLeft) < std::abs (rRight))
		m_MinVector = m_LeftBracket;
	else if (std::abs (rRight) < m_rTolerance)
		m_MinVector = m_RightBracket;
	else
		m_MinVector = vMiddleBracket;

	m_OOFValue = m_ptrOptimisableFunction->Evaluate(m_MinVector);

	return;
}

//////////////////////////////////////////////////
