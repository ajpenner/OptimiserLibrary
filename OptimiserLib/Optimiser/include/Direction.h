// Direction.h
//
// 2016
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#ifdef _WIN32
#include "../OptimiserDLL/Exporter.h"
#else
#define OPTIMISER
#endif
#include <armadillo>

class IFunction;

//////////////////////////////////////////////////

class OPTIMISER IDirection
{
public:
	virtual arma::vec GetDirection (const IFunction& function, const arma::vec& vPosition) const = 0;
	virtual ~IDirection () {};
};

//////////////////////////////////////////////////
