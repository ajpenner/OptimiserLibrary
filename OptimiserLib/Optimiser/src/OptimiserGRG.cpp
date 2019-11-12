// OptimiserGRG.cpp
//
// 2016
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif

//////////////////////////////////////////////////

#include "OptimiserGRG.h"
#include "optimisable.h"
#include "STLExtension.h"

// to find min/max of constraints
#include "MinimiserBrent.h"
#include "MaximiserBrent.h" // this should be a single class

#include <cassert>
#include <armadillo>

#include <boost/foreach.hpp>

//////////////////////////////////////////////////

void COptimiserGRG::AddConstraint (std::shared_ptr<const IConstraint> ptrConstraint)
{
	m_vConstraints.push_back (ptrConstraint);
}

//////////////////////////////////////////////////

void COptimiserGRG::SetUnconstrainedOptimiser (IOptimiser* pcOptimiser)
{
}

//////////////////////////////////////////////////

void COptimiserGRG::_CheckForInequalityConstraintsAndAddSlackVariables ()
{
	for (auto it = m_vConstraints.begin (); it != m_vConstraints.end (); ++it)
	{
		auto eType = (*it)->GetType();
		if ( eType == eUnknown)
		{
			// should throw
			return;
		}
		if ( eType == eEQ)
		{
			continue;
		}
		else
		{
			// added variable
			++m_iSlackCount;
		}
	}
}

//////////////////////////////////////////////////

void COptimiserGRG::SetBounds(std::vector<arma::vec>& vBounds)
{
	m_vCoordinateBounds.resize (vBounds.size ());
	std::copy (vBounds.cbegin (), vBounds.cend (), m_vCoordinateBounds.begin ());
}

//////////////////////////////////////////////////

arma::vec COptimiserGRG::FindFeasible () const
{
	int i = 0;
	auto vSlack = arma::vec (m_iSlackCount).zeros();
	for (auto it = m_vConstraints.begin (); it != m_vConstraints.end (); ++it, ++i) // this should be zip iterated
	{
		auto eType = (*it)->GetType();
		if (eType == eEQ)
			--i; // equality constraints do not need slack variables

		auto rValue = (*it)->Evaluate(m_MinVector);

		if ( eType == eLT || eType == eLE )
		{
			vSlack[i] = -rValue;
		}
		else
		{
			vSlack[i] = rValue; // slack variables are always positive
		}
	}

	return std::move(vSlack);
}

//////////////////////////////////////////////////

void COptimiserGRG::Optimise (double target)
{
	assert (m_ptrOptimisableFunction);

	// add slack variables
	_CheckForInequalityConstraintsAndAddSlackVariables();

	const auto iOriginalDim = m_ptrOptimisableFunction->GetFunctionDimension();

	auto iDim = iOriginalDim + m_iSlackCount;
	arma::vec yQuality(iDim); // rename

	auto vFeasibleSlackValues = FindFeasible();

	// merge original guess with slack variables
	m_MinVector.resize(iDim);
	for (auto i = iOriginalDim; i < iDim; ++i)
	{
		m_MinVector[i] = vFeasibleSlackValues (i - iOriginalDim);
	}

	// set bounds with slack variables
	m_vCoordinateMin.resize (iDim);
	m_vCoordinateMin.zeros ();
	m_vCoordinateMax.resize (iDim);
	m_vCoordinateMax.zeros ();

	auto i = arma::uword(0);
	for (; i < static_cast<arma::uword>(m_vCoordinateBounds.size()); ++i) // puke
	{
		m_vCoordinateMin[i] = m_vCoordinateBounds[i][0];
		m_vCoordinateMax[i] = m_vCoordinateBounds[i][1];
	}

	// use brent method to find min for the constraint functions inside the domain specified by the user
	CMinimiserBrent optimiser;
	optimiser.SetLeftBracket (m_vCoordinateMin); // should be 0,0
	optimiser.SetRightBracket (m_vCoordinateMax); // should be 5,5
	optimiser.SetTolerance (m_rTolerance);

	// use brent method to find max for the constraint functions inside the domain specified by the user
	CMaximiserBrent optimiser2;
	optimiser2.SetLeftBracket (m_vCoordinateMin); // should be 0,0
	optimiser2.SetRightBracket (m_vCoordinateMax); // should be 5,5
	optimiser2.SetTolerance (m_rTolerance);

	for (auto it = m_vConstraints.cbegin (); it != m_vConstraints.cend (); ++it, ++i)
	{
//		optimiser.SetOptimisableFunction(*it); // needs to be a shared_ptr
		optimiser.Optimise (0);
		auto minVec = optimiser.GetVector ();
		m_vCoordinateMin[i] = optimiser.GetOOF ();

//		optimiser2.SetOptimisableFunction (*it); // needs to be a shared_ptr
		optimiser2.Optimise (0);
		auto maxVec = optimiser2.GetVector ();
		m_vCoordinateMax[i] = optimiser2.GetOOF ();
	}

	for (auto it = m_vConstraints.cbegin (); it != m_vConstraints.cend (); ++it, ++i)
	{
		m_vCoordinateMin[i] = (*it)->Evaluate (m_vCoordinateMin);
		m_vCoordinateMax[i] = (*it)->Evaluate (m_vCoordinateMax);
	}

	// choose the slack variables
	ZippedTransform(m_MinVector, m_vCoordinateMin, m_vCoordinateMax, yQuality, 
		[] (const decltype(*m_MinVector.cbegin ())& rValue, const decltype(*m_vCoordinateMin.cbegin ())& rLower, const decltype(*m_vCoordinateMax.cbegin ())& rUpper)
	{ 
		double rTolerance = 1e-12;
		assert (rUpper - rLower > rTolerance);
		return std::min ((rValue-rLower), (rUpper-rValue)) / (rUpper-rLower); 
	});

	const auto cPermutation = sort_permutation<arma::uword>(yQuality, std::greater_equal<double>());
	yQuality = apply_permutation<arma::uword>(yQuality, cPermutation);

	// need to sort the various vectors according to the permutation, may want to just apply teh permutation in place each time. Will need to sort iteratively and this will be slow

	auto iK = static_cast<arma::uword>(m_vConstraints.size());
	auto vBasicVariables = yQuality.submat(1,1,1,iK);
	auto vNonBasicVariables = yQuality.submat(1,iK,1,iDim);

	// calculate the derivatives of the objective function these are also sorted by yQuality



	// calculate derivatives of the constraint equations with the original arrangement of the vectors

	// need a template specialisation overload to handle crowbegin etc...
//	m_mDerivatives = apply_permutation<arma::uword>( m_mDerivatives, cPermutation); // permute to match the yQuality arrangement

	auto m_mJDerivative = arma::mat (iK, iK); // == dh_(1->k) / dx_(1->k)
	
	// need to use matrix view to do this, or it will become too expensive
	// copy iK x iK elements from m_mDerivative to m_mJDer... 
	auto m_mCDerivative = arma::mat (iK, iDim-iK); // == dh_(1->k) / dx_(k->iDim)
	// copy iK x (iDim - iK) elements from m_mDerivative to m_mCDer... 




	return;
}

//////////////////////////////////////////////////
