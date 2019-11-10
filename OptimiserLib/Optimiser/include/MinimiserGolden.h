// MinimiserGolden.h
//
// 2016
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "optimiser.h"
#include "optimisable.h"
#include <limits>
#include <armadillo>

//////////////////////////////////////////////////

class OPTIMISER CMinimiserGolden : public BExtremum
{
public:
	CMinimiserGolden () {}
	~CMinimiserGolden () {}

	void SetLeftBracket (const arma::vec& v) override { m_LeftBracket = v; }
	void SetRightBracket (const arma::vec& v) override { m_RightBracket = v; }

	void SetDirectionVector (arma::vec& v) override { m_MinVector = v; }
	void SetLeftBracketRange (double rLeft) override { m_LeftBracket = (rLeft*m_MinVector).eval (); }
	void SetRightBracketRange (double rRight) override { m_RightBracket = (rRight*m_MinVector).eval (); }

	void SetParameters (const std::vector<variant>& parameters) override {};

	void Optimise (double /*target*/) override;

private:

	class Transform : public IFunction
	{
	private:
		std::shared_ptr<const IFunction> m_ptrMainFunction;
		double m_lBracket;
		double m_rBracket;

	public:
		Transform (std::shared_ptr<const IFunction> mainFunction, double lBracket, double rBracket) : 
			m_ptrMainFunction (mainFunction), 
			m_lBracket(lBracket), 
			m_rBracket(rBracket)
		{
		}

		double Evaluate (const arma::vec& vValues) const override // input is w
		{
			arma::vec x = vValues*(m_rBracket - m_lBracket) + m_lBracket; // function evaluates at x
			return m_ptrMainFunction->Evaluate (x);
		}
		
		arma::vec CalculateGradient (const arma::vec& vec) const
		{
			arma::vec x = vec*(m_rBracket - m_lBracket) + m_lBracket; // function evaluates at x
			return m_ptrMainFunction->CalculateGradient (x);
		};

		arma::mat CalculateHessian (const arma::vec& vec) const
		{
			arma::vec x = vec*(m_rBracket - m_lBracket) + m_lBracket; // function evaluates at x
			return m_ptrMainFunction->CalculateHessian (x);
		};
		
		arma::uword GetFunctionDimension () const
		{
			return m_ptrMainFunction->GetFunctionDimension ();
		};

		arma::vec GetOriginalCoordinate (double w) const
		{
			return arma::vec ( {w*(m_rBracket - m_lBracket) + m_lBracket} );
		};
	};

	void CalculateUpperBracket (arma::vec& vUpperBracket);
	void CalculateLowerBracket (arma::vec& vLowerBracket);

	arma::vec m_LeftBracket;
	arma::vec m_RightBracket;

	const double m_rGoldenRatio = 1.0 - (3.0-std::sqrt (5.0)) / 2.0;
};

//////////////////////////////////////////////////