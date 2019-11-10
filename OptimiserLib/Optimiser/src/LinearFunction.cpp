// VectorFunction.cpp
//
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "stdafx.h"
#include "LinearFunction.h"

//////////////////////////////////////////////////

CLinearFunction::CLinearFunction (std::shared_ptr<const IFunction> mainFunction) : m_ptrMainFunction(mainFunction)
{
}

//////////////////////////////////////////////////

void CLinearFunction::SetVector (const arma::vec& vector)
{
	m_Vector = vector;
	m_iDimension = static_cast<uint8_t>(m_Vector.size());
}

//////////////////////////////////////////////////

void CLinearFunction::SetOffset (const arma::vec& vector)
{
	m_Offset = vector;
}

//////////////////////////////////////////////////

double CLinearFunction::Evaluate(const arma::vec& vValues) const
{
	assert (vValues.n_cols == 1);
	arma::vec X = m_Offset + vValues.at(0)*m_Vector;
	return m_ptrMainFunction->Evaluate(X);
}

//////////////////////////////////////////////////

arma::vec CLinearFunction::CalculateGradient (const arma::vec& vec) const
{ 
	return arma::vec (m_iDimension).zeros (); 
}

//////////////////////////////////////////////////

arma::mat CLinearFunction::CalculateHessian (const arma::vec& vec) const
{ 
	return arma::mat (m_iDimension, m_iDimension).zeros (); 
}

//////////////////////////////////////////////////

arma::uword CLinearFunction::GetFunctionDimension () const
{ 
	return m_iDimension; 
}

//////////////////////////////////////////////////