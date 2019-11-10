// OptimiserComplex.cpp
//
// 2017
//////////////////////////////////////////////////

#include "stdafx.h"

//////////////////////////////////////////////////

#include "OptimiserComplex.h"
#include "optimisable.h"
#include "CompoundMOM.h"

#include <cassert>
#include <iostream>
#include <armadillo>

#include <boost/foreach.hpp>
#include "boost\iterator\zip_iterator.hpp"


//////////////////////////////////////////////////

void COptimiserComplex::AddConstraint ( std::shared_ptr<const IConstraint> ptrConstraint )
{
	m_vConstraints.push_back(ptrConstraint);
}

//////////////////////////////////////////////////

void COptimiserComplex::SetParameters (const std::vector<variant>& parameters)
{
	m_beta = boost::get<double> (parameters[eComplexParameters::eBeta]);
	assert (m_beta >= 0.0 && m_beta <= 1.0);

	m_gamma = boost::get<double> (parameters[eComplexParameters::eGamma]);
	assert (m_gamma >= 1.0);

	assert (m_ptrOptimisableFunction);
	for (int i = 0; i < m_ptrOptimisableFunction->GetFunctionDimension (); ++i)
	{
		auto pr = boost::get<std::pair<double,double>> (parameters[eComplexParameters::eRange + i]);
		m_vRange.push_back(pr);
	}
}

//////////////////////////////////////////////////

bool COptimiserComplex::VerifyConstraints(const arma::vec& vec)
{
	auto satisfied = false;
	BOOST_FOREACH (const auto& constraint, m_vConstraints)
	{
		satisfied = constraint->IsSatisfied (vec);
		if (!satisfied) break;
	}
	return satisfied;
}

//////////////////////////////////////////////////

void COptimiserComplex::EnforceRange(arma::vec& vec)
{
	assert (vec.size () == m_dim);
	for(uint8_t i =0; i < m_dim; ++i)
	{
		auto& value = vec[i];
		if (value < m_vRange[i].first)
			value = m_vRange[i].first;
		else if (value > m_vRange[i].second)
			value = m_vRange[i].second;
	}
}

//////////////////////////////////////////////////

void COptimiserComplex::EnforceConstraints(arma::vec& vec)
{
	assert (vec.size () == m_dim);

}

//////////////////////////////////////////////////

void COptimiserComplex::SetInitialPoints()
{
	auto& distribution = *m_ptrDistribution;
	for (auto i = 0; i < m_simplexCount; ++i)
	{
		auto satisfied = false;
		do
		{
			arma::vec coord = arma::zeros(m_dim);
			for (size_t i = 0; i < m_dim; ++i)
			{
				auto range = m_vRange.at(i);
				auto x = range.first + distribution (*m_ptrGenerator)*(range.second - range.first);
				coord[i] = x;
			}
			satisfied = VerifyConstraints (coord);
			if (satisfied) m_vCoordinates[i] = coord;

		} while (!satisfied);
	}
}

//////////////////////////////////////////////////

void COptimiserComplex::Optimise (double target)
{
	assert (m_ptrOptimisableFunction);
	assert (m_ptrUnconstrainedOptimiser);

	m_dim = m_ptrOptimisableFunction->GetFunctionDimension();
	m_simplexCount = 2 * m_ptrOptimisableFunction->GetFunctionDimension ();
	// create 2N initial feasible points
	m_vCoordinates.resize (m_simplexCount);
	m_evaluations.resize (m_simplexCount);

	assert (m_vRange.size () == m_ptrOptimisableFunction->GetFunctionDimension ());
	SetInitialPoints ();

	std::transform (m_vCoordinates.cbegin (), m_vCoordinates.cend (), m_evaluations.begin (),
		[&](const auto& v)
	{
		return m_ptrOptimisableFunction->Evaluate (v);
	});

	auto escape = realEmpty;
	do
	{
		auto itWorst = std::max_element(m_evaluations.cbegin(), m_evaluations.cend());
		arma::vec xW = m_vCoordinates[std::distance(m_evaluations.cbegin(), itWorst)];
		auto fW = *itWorst;
		m_vCoordinates.erase(m_vCoordinates.begin() + std::distance(m_evaluations.cbegin(), itWorst));
		m_evaluations.erase(m_evaluations.begin() + std::distance(m_evaluations.cbegin(), itWorst));

		auto itWorstRemaining = std::max_element(m_evaluations.begin(), m_evaluations.end());
//		auto thing = std::distance(m_evaluations.begin(), itWorstRemaining);
		auto& xG = m_vCoordinates[std::distance(m_evaluations.begin(), itWorstRemaining)];
		auto fG = *itWorstRemaining;

		auto itBest = std::min_element(m_evaluations.cbegin(), m_evaluations.cend());
		auto& xB = m_vCoordinates[std::distance(m_evaluations.cbegin(), itBest)];
		auto fB = *itBest;

		// obtain the centroid
		arma::vec zero = arma::zeros (m_dim);
		arma::vec centroid = std::accumulate (m_vCoordinates.cbegin (), m_vCoordinates.cend (), zero,
			[](const auto& x, const auto& y) {return x + y; }) / static_cast<double>(m_dim);

		auto fC = m_ptrOptimisableFunction->Evaluate(centroid);

		// reflect
		arma::vec xR = 2*centroid - xW;
		EnforceRange (xR);

		auto xNew = xR;

		auto fR = m_ptrOptimisableFunction->Evaluate (xR);
		if (fR < fB) // stretch
		{
			xNew = (1 + m_gamma)*centroid - m_gamma*xW;
		}
		else if (fR >= fW) // retract
		{
			xNew = (1 - m_beta)*centroid + m_beta*xW;
		}
		else if (fG < fR && fR < fW) // retract
		{
			xNew = (1 + m_beta)*centroid - m_beta*xW;
		}
		EnforceRange (xNew);
		EnforceConstraints(xNew);
		// need to know that new point is feasible
		bool thing = VerifyConstraints (xNew);

		while (!VerifyConstraints(xNew))
		{
			xNew = (1 - m_beta)*centroid + m_beta*xNew;
		}

		m_vCoordinates.push_back (xNew);
		m_evaluations.push_back (m_ptrOptimisableFunction->Evaluate (xNew));

		// check for termination
		escape = std::accumulate (m_evaluations.cbegin (), m_evaluations.cend (), 0.0,
			[](auto x, auto y) { return x + std::pow(y,2); }) / static_cast<double>(m_dim);
		// also need to check delta on the coordinates
	} while (escape > m_rTolerance);

	return;
}

//////////////////////////////////////////////////