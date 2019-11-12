// DirectionQuasiNewton.cpp
//
// 2016
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif

//////////////////////////////////////////////////

#include "BackTracking.h"
#include "optimisable.h"
#include "DirectionQuasiNewton.h"

//////////////////////////////////////////////////

arma::vec CDirectionQuasiNewton::GetDirection (const IFunction& function, const arma::vec& vPosition) const
{
	auto vGrad = function.CalculateGradient (vPosition);
	// -H*vGrad where H is an approximation to an inverse Hessian
	throw std::runtime_error("Not yet implemented");
	return vGrad;// (-arma::inv_sympd (mHess)*vGrad).eval;
}

//////////////////////////////////////////////////
