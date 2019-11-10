// OptimiserEllipsoid.cpp
//
// 2017
//////////////////////////////////////////////////

#include "stdafx.h"

//////////////////////////////////////////////////

#include "OptimiserEllipsoid.h"
#include "optimisable.h"
#include <cassert>
#include <iostream>
#include <armadillo>

#include <boost/foreach.hpp>

//////////////////////////////////////////////////

void COptimiserEllipsoid::AddConstraint ( std::shared_ptr<const IConstraint> ptrConstraint )
{
	auto ptrLinearConstraint = std::dynamic_pointer_cast<const ILinearConstraint>(ptrConstraint);
	m_vConstraints.push_back(ptrLinearConstraint);
}

//////////////////////////////////////////////////

void COptimiserEllipsoid::SetUnconstrainedOptimiser (IOptimiser* pcOptimiser)
{
}

//////////////////////////////////////////////////

void COptimiserEllipsoid::SetParameters (const std::vector<variant>& parameters)
{
	m_rRadius = boost::get<double> (parameters[eEllipsoidParameters::eRadius]);
	assert (m_rRadius >= realZero);
}

//////////////////////////////////////////////////

void COptimiserEllipsoid::Optimise (double target)
{
	assert (m_ptrOptimisableFunction);

	const auto iDim = m_ptrOptimisableFunction->GetFunctionDimension ();
	arma::vec a(iDim); // rename
	auto Q = arma::mat(iDim, iDim).eye();
	auto C = arma::vec(iDim).zeros(); // rename

	// set radius
	Q *= m_rRadius*m_rRadius;

	auto rDim = static_cast<double>(iDim);
	size_t uIteration = 0;
	while( m_uMaxIterati Sons > uIteration)
	{
		auto violatedConstraint = std::find_if (m_vConstraints.cbegin (), m_vConstraints.cend (), 
			[this](const decltype(*m_vConstraints.cbegin ())& value) 
		{ 
			return !value->IsSatisfied (this->GetVector ()); 
		});

		if (violatedConstraint == m_vConstraints.cend())
		{
			return;
		}

		a = (*violatedConstraint)->GetCoeff ();

		auto v = ((Q*a) / sqrt ((a.t ()*Q*a).eval ()[0])).eval ();

		m_MinVector -= v / (rDim + 1);

		Q = (rDim*rDim) / (rDim*rDim - 1)*(Q - 2 / (rDim + 1)*v*v.t ()).eval ();

		++uIteration;
	}

	return;
}

//////////////////////////////////////////////////