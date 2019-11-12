// WeakWolfe.h
//
// 2016
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "LineSearch.h"
#include <armadillo>

//////////////////////////////////////////////////

class IFunction;
class IDirection;

//////////////////////////////////////////////////

class OPTIMISER CLineSearchWeakWolfe : public ILineSearch
{
public:
	CLineSearchWeakWolfe () 
		: m_rTolerance(std::numeric_limits<double>::max())
		, m_rStepLength(1.0)
		, m_ptrFunction(nullptr) {}

	~CLineSearchWeakWolfe () {}

	void SetFunction (std::shared_ptr<const IFunction> ptrFunction) override;
	void SetDirectionMethod (const IDirection* pcDirection) override;
	void SetTolerance (double rTolerance) override;
	void CalculateStepLength (const arma::vec& vPosition) const override;

	double GetStepLength () const override;

private:

	double m_rTolerance;
	mutable double m_rStepLength;
	std::shared_ptr<const IFunction> m_ptrFunction;
	const IDirection* m_pcDirection;
};

//////////////////////////////////////////////////
