// StrongWolfe.cpp
//
// 2016
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif

//////////////////////////////////////////////////

#include "StrongWolfe.h"
#include "optimisable.h"
#include "Direction.h"

//////////////////////////////////////////////////

void CLineSearchStrongWolfe::SetFunction (std::shared_ptr<const IFunction> ptrFunction)
{
	m_ptrFunction = ptrFunction;
}

//////////////////////////////////////////////////

void CLineSearchStrongWolfe::SetDirectionMethod (const IDirection* pcDirection)
{
	m_pcDirection = pcDirection;
}

//////////////////////////////////////////////////

void CLineSearchStrongWolfe::SetTolerance (double rTolerance)
{
	m_rTolerance = rTolerance;
}

//////////////////////////////////////////////////

double CLineSearchStrongWolfe::Zoom (double rAlphaLo, double rAlphaHi,
	const arma::vec& vPosition,
	const arma::vec& vDirection,
	double rFunction,
	double rFp) const
{
	auto rAlpha = 0.0;
	do
	{
		// interpolate between the input arguments;
		// for now I use the midpoint
		rAlpha = 0.5*(rAlphaHi + rAlphaLo);
		arma::vec vPositionp1 = (vPosition + rAlpha*vDirection).eval ();
		auto rFunctionp1 = m_ptrFunction->Evaluate (vPositionp1);

		arma::vec vPositionLo = (vPosition + rAlphaLo*vDirection).eval ();
		auto rFunctionLo = m_ptrFunction->Evaluate (vPositionLo);

		if ((rFunctionp1 > rFunction + m_rC1*rAlpha*rFp)
			|| (rFunction >= rFunctionLo))
		{
			rAlphaHi = rAlpha; // updated local hi
		}
		else
		{
			auto vGradp1 = m_ptrFunction->CalculateGradient (vPositionp1);
			auto rFp1p = arma::dot (vGradp1, vDirection);

			if (std::abs (rFp1p) <= -m_rC2*rFp)
			{
				break;
			}

			if (rFp1p*(rAlphaHi - rAlphaLo) >= 0.0)
			{
				rAlphaHi = rAlphaLo;
			}
			rAlphaLo = rAlpha;
		}
		if(rAlphaHi == rAlphaLo)
		{
			break;
		}
	} while (true);

	assert (rAlpha > 0.0);
	return rAlpha;
}

//////////////////////////////////////////////////

void CLineSearchStrongWolfe::CalculateStepLength (const arma::vec& vPosition) const
{
	auto rAlphaOld = 0.0;
	auto rAlphaMax = 1; // is 1 sensible?
	m_rStepLength = 0.4; //(between 0 and rAlphaMax), maybe random?

	auto rFunction = m_ptrFunction->Evaluate (vPosition); // phi(0)
	auto vGrad = m_ptrFunction->CalculateGradient (vPosition);
	const auto vDirection = m_pcDirection->GetDirection (*m_ptrFunction, vPosition);
	auto rFp = arma::dot (vGrad, vDirection); // phi'(0)

	auto rFunctionOld = rFunction;

	auto iCount = 0;
	do
	{
		arma::vec vPositionp1 = (vPosition + m_rStepLength*vDirection).eval();

		auto rFunctionp1 = m_ptrFunction->Evaluate (vPositionp1);

		if ( (rFunctionp1 > rFunction + m_rC1*m_rStepLength*rFp)
			|| (iCount > 0 && rFunctionp1 >= rFunctionOld) )
		{
			m_rStepLength = Zoom (rAlphaOld, m_rStepLength, vPosition, vDirection, rFunction, rFp);
			break;
		}

		auto vGradp1 = m_ptrFunction->CalculateGradient (vPositionp1);
		auto rFp1p = arma::dot (vGradp1, vDirection);

		if (std::abs(rFp1p) < -m_rC2*rFp)
		{
			break;
		}

		if(rFp1p >= 0.0)
		{
			m_rStepLength = Zoom (m_rStepLength, rAlphaOld, vPosition, vDirection, rFunction, rFp);
			break;
		}

		rAlphaOld = m_rStepLength;
		m_rStepLength += 0.01; // between (m_rStepLength, rAlphaMax)
		++iCount;
	} while (true); // forever
}

//////////////////////////////////////////////////

double CLineSearchStrongWolfe::GetStepLength () const
{
	return m_rStepLength;
}

//////////////////////////////////////////////////
