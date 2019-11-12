// DirectionCauchy.cpp
//
// 2016
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif

//////////////////////////////////////////////////

#include "BackTracking.h"
#include "optimisable.h"
#include "DirectionCauchy.h"

//////////////////////////////////////////////////

arma::vec CDirectionCauchy::GetDirection(const IFunction& function, const arma::vec& vPosition) const
{
	return -arma::normalise(function.CalculateGradient(vPosition));
}

//////////////////////////////////////////////////
