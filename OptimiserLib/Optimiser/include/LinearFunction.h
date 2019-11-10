// VectorFunction.h
//
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "optimisable.h"
#include <memory>

//////////////////////////////////////////////////

class FUNCTION CLinearFunction : public IFunction
{
private:
	uint8_t m_iDimension;
	arma::vec m_Vector; // needs DLL interface (but never actually used across a DLL)
	arma::vec m_Offset;
	std::shared_ptr<const IFunction> m_ptrMainFunction;

public:

	CLinearFunction (std::shared_ptr<const IFunction> mainFunction);
	~CLinearFunction() {}

	void SetVector (const arma::vec& vector);
	void SetOffset (const arma::vec& vector);

	double Evaluate (const arma::vec& vValues) const override;

	arma::vec CalculateGradient (const arma::vec& vec) const override;
	arma::mat CalculateHessian (const arma::vec& vec) const override;

	arma::uword GetFunctionDimension () const;
};

//////////////////////////////////////////////////