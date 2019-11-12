// StrongWolfe.h
//
// 2016
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "LineSearch.h"
#include <armadillo>
#include <cassert>

//////////////////////////////////////////////////

class IFunction;
class IDirection;

//////////////////////////////////////////////////

class OPTIMISER CLineSearchStrongWolfe : public ILineSearch
{
public:
	CLineSearchStrongWolfe () 
		: m_rTolerance (std::numeric_limits<double>::max())
		, m_rStepLength (1.0)
		, m_ptrFunction (nullptr) 
	{
		m_rC1 = 0.01; // set by user...
		assert (0.0 < m_rC1);

		m_rC2 = 0.5;
		assert (m_rC2 < 1.0);
		assert (m_rC1 < m_rC2);
	};
	~CLineSearchStrongWolfe () {};

	void SetFunction (std::shared_ptr<const IFunction> ptrFunction) override;
	void SetDirectionMethod (const IDirection* pcDirection) override;
	void SetTolerance (double rTolerance) override;
	void CalculateStepLength (const arma::vec& vPosition) const override;

	double GetStepLength () const override;

private:

	double Zoom (double rAlphaLo, double rAlphaHi,
		const arma::vec& vPosition,
		const arma::vec& vDirection,
		double rFunction,
		double rFp) const;

	double m_rTolerance;
	mutable double m_rStepLength;
	std::shared_ptr<const IFunction> m_ptrFunction;
	const IDirection* m_pcDirection;

	double m_rC1;
	double m_rC2;
};

//////////////////////////////////////////////////
