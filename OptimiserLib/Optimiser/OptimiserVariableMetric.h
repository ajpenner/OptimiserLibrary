// OptimiserVariableMetric.h
//
// 2016
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "optimiser.h"
#include "LineSearch.h"

//////////////////////////////////////////////////

class OPTIMISER IOptimiserVariableMetric : public BOptimiser
{
public:

	virtual void SetLineSearchMethod( std::shared_ptr<ILineSearch> ptrLineSearch ) = 0;
};

//////////////////////////////////////////////////

class OPTIMISER BOptimiserVariableMetric : public IOptimiserVariableMetric
{
public:

	void SetLineSearchMethod (std::shared_ptr<ILineSearch> ptrLineSearch )
	{
		m_ptrLineSearch = ptrLineSearch;
	}

protected: // function
	virtual arma::mat UpdateHessian() = 0;

protected: // members
	std::shared_ptr<ILineSearch> m_ptrLineSearch;

	arma::mat m_mEye;
	arma::mat m_mHess;
	arma::vec m_deltaX;
	arma::vec m_deltaGrad;
};

//////////////////////////////////////////////////