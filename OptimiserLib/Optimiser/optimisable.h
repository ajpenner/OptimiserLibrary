// optimisable.h
//
// 2016
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "../FunctionDLL/Exporter.h"
#include <armadillo>

//////////////////////////////////////////////////

class FUNCTION IFunction // need to default numeric derivatives, override can give us closed-form 
{
public:
	virtual double Evaluate( const arma::vec& vec ) const = 0;
	virtual arma::vec CalculateGradient( const arma::vec& vec ) const = 0;
	virtual arma::mat CalculateHessian( const arma::vec& vec ) const = 0;
	virtual arma::uword GetFunctionDimension() const = 0;

	virtual ~IFunction () {};
};

//////////////////////////////////////////////////