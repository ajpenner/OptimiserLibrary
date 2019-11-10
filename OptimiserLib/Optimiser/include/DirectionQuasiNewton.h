// DirectionQuasiNewton.h
//
// 2016
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include <armadillo>
#include "Direction.h"
#include "optimisable.h"

//////////////////////////////////////////////////

class OPTIMISER CDirectionQuasiNewton : public IDirection
{
public:
	arma::vec GetDirection (const IFunction& function, const arma::vec& vPosition) const override;
};

//////////////////////////////////////////////////