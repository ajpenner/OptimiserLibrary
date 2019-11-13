// constraintsImpl.h
//
// 2016
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif

//////////////////////////////////////////////////

#include "constraints.h"
#include <armadillo>
#include <functional>

//////////////////////////////////////////////////

class CLinearConstraintImpl : public ILinearConstraint
{
private:

	arma::vec m_vConstraintCoeff;
	double m_rConstraintLimit;
	EConstraintType m_eType;

public:

	CLinearConstraintImpl (EConstraintType ePredicate, double rLimit, const arma::vec& vCoeff) 
		: m_vConstraintCoeff (vCoeff)
		, m_rConstraintLimit (rLimit)
		, m_eType (ePredicate) 
	{
	}

	arma::vec GetCoeff () const override { return m_vConstraintCoeff; }
	void SetCoeff (arma::vec& vCoef) override { m_vConstraintCoeff = vCoef; }
	double GetLimit () const override { return m_rConstraintLimit; }
	void SetLimit (double rValue) override { m_rConstraintLimit = rValue; }
	double Evaluate (const arma::vec& vCoordinates) const override
	{
		return arma::dot(vCoordinates,m_vConstraintCoeff);
	}

	bool IsSatisfied (const arma::vec& vSolution) const 
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
			return Evaluate (vSolution) ==  m_rConstraintLimit;
			break;
		default:
			return false; // no predicate set, should throw
		}
	}
	// From IOptimisable, but not needed for this type of constraint (would rather throw)
	arma::vec CalculateGradient (const arma::vec& vec) const override { return arma::vec (GetFunctionDimension ()).zeros(); };
	arma::mat CalculateHessian (const arma::vec& vec) const override { return arma::mat (GetFunctionDimension (), GetFunctionDimension ()).zeros (); };
	arma::uword GetFunctionDimension () const override { return m_vConstraintCoeff.size (); };

	EConstraintType GetType () const override { return m_eType; }
};
