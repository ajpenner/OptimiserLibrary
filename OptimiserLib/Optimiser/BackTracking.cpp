// BackTracking.cpp
//
// 2016
//////////////////////////////////////////////////

#include "stdafx.h"

//////////////////////////////////////////////////

#include "BackTracking.h"
#include "optimisable.h"
#include "Direction.h"

//////////////////////////////////////////////////

void CLineSearchBackTrack::SetFunction(std::shared_ptr<const IFunction> ptrFunction)
{
	m_ptrFunction = ptrFunction; 
}

//////////////////////////////////////////////////

void CLineSearchBackTrack::SetDirectionMethod (const IDirection* pcDirection)
{
	m_pcDirection = pcDirection;
}

//////////////////////////////////////////////////

void CLineSearchBackTrack::SetTolerance (double rTolerance)
{
	m_rTolerance = rTolerance;
}

//////////////////////////////////////////////////
// this is not right
double CLineSearchBackTrack::_DFunction(const arma::vec& vPosition, const arma::vec& vDirection) const
{
	auto rFp1 = m_ptrFunction->Evaluate(vPosition + 0.01*vDirection);
	auto rF = m_ptrFunction->Evaluate (vPosition);

	return (rFp1 - rF) / 0.01;
}

//////////////////////////////////////////////////

void CLineSearchBackTrack::CalculateStepLength(const arma::vec& vPosition) const
{
	auto rGamma = 0.5e0; // [0.5, 0.8] // set by user
	auto rC = 1.0e-3; // [1e-3,1e-1]

	auto rFunction = m_ptrFunction->Evaluate(vPosition);
	auto vGrad = m_ptrFunction->CalculateGradient(vPosition);

	const auto vDirection = m_pcDirection->GetDirection(*m_ptrFunction, vPosition);
	auto rDeltaFunction = rC*arma::dot(vGrad,vDirection);

	auto rFunctionp1 = m_ptrFunction->Evaluate (vPosition+vDirection);

	do
	{
		m_rStepLength = rGamma*m_rStepLength;
		rFunctionp1 = m_ptrFunction->Evaluate(vPosition + m_rStepLength*vDirection);
	} while (rFunctionp1 > rFunction + m_rStepLength*rDeltaFunction);
}

//////////////////////////////////////////////////

double CLineSearchBackTrack::GetStepLength() const
{ 
	return m_rStepLength; 
}

//////////////////////////////////////////////////