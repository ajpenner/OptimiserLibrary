// CNonlinearConstraintImpl.h
//
// 2016
//////////////////////////////////////////////////

#include "stdafx.h"

//////////////////////////////////////////////////

#include "constraints.h"
#include <armadillo>
#include <memory>

//////////////////////////////////////////////////

class CNonlinearConstraintImpl : public IConstraint
{
private:

	// the RHS
	double m_rConstraintLimit; // for convenience default to 0
							   // compare LHS to RHS
	// the LHS
	std::shared_ptr<const IFunction> m_ptrConstraint;

	EConstraintType m_eType;

public:

	CNonlinearConstraintImpl(EConstraintType ePredicate, std::shared_ptr<const IFunction> ptrConstraint, double rLimit = realZero) :
		m_rConstraintLimit (rLimit), m_ptrConstraint(ptrConstraint), m_eType(ePredicate)
	{
	};

	~CNonlinearConstraintImpl () {};

	double GetLimit () const override { return m_rConstraintLimit; };
	void SetLimit (double rValue) override { m_rConstraintLimit = rValue; };
	bool IsSatisfied (const arma::vec& vSolution) const override
	{
		switch (m_eType)
		{
		case eGT:
			return Evaluate (vSolution) > m_rConstraintLimit;
			break;
		case eGE:
			return Evaluate (vSolution) >= m_rConstraintLimit;
			break;
		case eLT:
			return Evaluate (vSolution) < m_rConstraintLimit;
			break;
		case eLE:
			return Evaluate (vSolution) <= m_rConstraintLimit;
			break;
		case eEQ:
			return Evaluate (vSolution) == m_rConstraintLimit;
			break;
		default:
			return false; // no predicate set, should throw
		}
	}

	double Evaluate (const arma::vec& vCoordinates) const override
	{ 
		return m_ptrConstraint->Evaluate(vCoordinates); 
	};

	arma::vec CalculateGradient (const arma::vec& vCoordinates) const override { return m_ptrConstraint->CalculateGradient (vCoordinates); };
	arma::mat CalculateHessian (const arma::vec& vCoordinates) const override { return m_ptrConstraint->CalculateHessian(vCoordinates); };
	arma::uword GetFunctionDimension () const override { return m_ptrConstraint->GetFunctionDimension(); };

	EConstraintType GetType () const override { return m_eType; }
};