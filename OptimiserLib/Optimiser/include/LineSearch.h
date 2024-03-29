// LineSearch.h
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
#include <memory>

//////////////////////////////////////////////////

class IFunction;
class IDirection;

//////////////////////////////////////////////////

class OPTIMISER ILineSearch
{
public:
	virtual void SetFunction (std::shared_ptr<const IFunction> ptrFunction) = 0;
	virtual void SetDirectionMethod (const IDirection* pcDirection) = 0;
	virtual void SetTolerance (double rTolerance) = 0;
	virtual void CalculateStepLength (const arma::vec& vPosition) const = 0;

	virtual double GetStepLength () const = 0;
};

//////////////////////////////////////////////////
