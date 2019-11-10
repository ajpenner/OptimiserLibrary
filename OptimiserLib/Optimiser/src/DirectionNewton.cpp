// DirectionNewton.cpp
//
// 2016
//////////////////////////////////////////////////

#include "stdafx.h"

//////////////////////////////////////////////////

#include "BackTracking.h"
#include "optimisable.h"
#include "DirectionNewton.h"

//////////////////////////////////////////////////

arma::vec CDirectionNewton::GetDirection(const IFunction& function, const arma::vec& vPosition) const
{
	auto mHess = function.CalculateHessian(vPosition);
	auto vGrad = function.CalculateGradient(vPosition);
	return (-arma::inv_sympd(mHess)*vGrad).eval();
}

//////////////////////////////////////////////////