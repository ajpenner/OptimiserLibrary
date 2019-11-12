// WeakWolfe.cpp
//
// 2016
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif

//////////////////////////////////////////////////

#include "WeakWolfe.h"
#include "optimisable.h"
#include "Direction.h"
#include <cassert>

//////////////////////////////////////////////////

void CLineSearchWeakWolfe::SetFunction (std::shared_ptr<const IFunction> ptrFunction)
{
	m_ptrFunction = ptrFunction;
}

//////////////////////////////////////////////////

void CLineSearchWeakWolfe::SetDirectionMethod (const IDirection* pcDirection)
{
	m_pcDirection = pcDirection;
}

//////////////////////////////////////////////////

void CLineSearchWeakWolfe::SetTolerance (double rTolerance)
{
	m_rTolerance = rTolerance;
}

//////////////////////////////////////////////////

void CLineSearchWeakWolfe::CalculateStepLength (const arma::vec& vPosition) const
{
	auto rC1 = 0.01;
	assert (0.0 < rC1 );

	auto rC2 = 0.5;
	assert (rC2 < 1.0);
	assert (rC1 < rC2);

	auto rAlpha = 0.0;
	m_rStepLength = 1.0;
	auto rBeta = std::numeric_limits<double>::max();

	auto rFunction = m_ptrFunction->Evaluate(vPosition);
	auto vGrad = m_ptrFunction->CalculateGradient(vPosition);
	const auto vDirection = m_pcDirection->GetDirection(*m_ptrFunction, vPosition);
	auto rFp = arma::dot (vGrad, vDirection);

	do
	{
		arma::vec vPositionp1 = (vPosition + m_rStepLength*vDirection).eval ();

		auto rFunctionp1 = m_ptrFunction->Evaluate (vPositionp1);

		auto vGradp1 = m_ptrFunction->CalculateGradient (vPositionp1);
		auto rFp1p = arma::dot (vGradp1, vDirection);

		if (rFunctionp1 > rFunction + rC1*m_rStepLength*rFp)
		{
			rBeta = m_rStepLength;
			m_rStepLength = 0.5*(rAlpha + rBeta);
		}
		else if (rFp1p < rC2*rFp)
		{
			rAlpha = m_rStepLength;
			if (rBeta == std::numeric_limits<double>::max())
			{
				m_rStepLength = 2.0*rAlpha;
			}
			else
			{
				m_rStepLength = 0.5*(rAlpha + rBeta);
			}
		}
		else
		{
			break;
		}

	} while (true); // forever
}

//////////////////////////////////////////////////

double CLineSearchWeakWolfe::GetStepLength () const
{
	return m_rStepLength;
}

//////////////////////////////////////////////////
