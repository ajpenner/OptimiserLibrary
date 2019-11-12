// OptimiserRandom.cpp
//
// 2017
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif

//////////////////////////////////////////////////

#include "OptimiserRandom.h"
#include "optimisable.h"
#include "CompoundMOM.h"

#include <cassert>
#include <armadillo>

#include <boost/foreach.hpp>
#include <boost/iterator/zip_iterator.hpp>


//////////////////////////////////////////////////

void COptimiserRandom::AddConstraint ( std::shared_ptr<const IConstraint> ptrConstraint )
{
	m_vConstraints.push_back(ptrConstraint);
}

//////////////////////////////////////////////////

void COptimiserRandom::SetParameters (const std::vector<variant>& parameters)
{
	m_rScale = boost::get<double> (parameters[eRandomParameters::eScale]);
	assert (m_rScale < 1.0);

	m_iPointCount = boost::get<int> (parameters[eRandomParameters::ePointCount]);
	assert (m_iPointCount > 0);

	assert (m_ptrOptimisableFunction);
	for (int i = 0; i < m_ptrOptimisableFunction->GetFunctionDimension (); ++i)
	{
		auto pr = boost::get<std::pair<double,double>> (parameters[eRandomParameters::eRange + i]);
		m_vRange.push_back(pr);
	}
}

//////////////////////////////////////////////////

bool COptimiserRandom::VerifyConstraints(const arma::vec& vec)
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

void COptimiserRandom::EnforceRange(arma::vec& vec)
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

void COptimiserRandom::SetRandomPoints()
{
	auto& distribution = *m_ptrDistribution;
	for (auto i = 0; i < m_iPointCount; ++i)
	{
		auto satisfied = false;
		do
		{ // I have hard coded 2 dimensions
			auto range = m_vRange.at (0); // this is terrible... I am tired
			auto x = m_MinVector.at(0) + range.first + distribution (*m_ptrGenerator)*(range.second-range.first);
			range = m_vRange.at (1);
			auto y = m_MinVector.at(1) + range.first + distribution (*m_ptrGenerator)*(range.second - range.first);
			auto coord = arma::vec ({ x,y });
			satisfied = VerifyConstraints (coord);
			if (satisfied) m_vCoordinates.push_back(coord);

		} while (!satisfied);
	}
}

//////////////////////////////////////////////////

void COptimiserRandom::Optimise (double target)
{
	assert (m_ptrOptimisableFunction);

	m_dim = m_ptrOptimisableFunction->GetFunctionDimension();

	assert (m_vRange.size () == m_dim);

	auto escape = std::numeric_limits<double>::max();
	do
	{
		SetRandomPoints ();
		m_evaluations.resize (m_vCoordinates.size());

		std::transform (m_vCoordinates.cbegin (), m_vCoordinates.cend (), m_evaluations.begin (),
			[&](const auto& v)
		{
			return m_ptrOptimisableFunction->Evaluate (v);
		});

		// find best point
		{
			auto best = m_evaluations.begin();
			for (auto it = m_evaluations.begin(); it != m_evaluations.end(); ++it)
			{
				if (*it < *best)
					best = it;
			}
			// remove the worst element
			auto distance = std::distance(m_evaluations.begin(), best);
			auto bestVector = *(m_vCoordinates.begin () + distance);
			auto bestOOF = *(m_evaluations.begin () + distance);
			m_vCoordinates.clear();
			m_evaluations.clear();

			if (bestOOF < m_OOFValue)
			{
				m_MinVector = bestVector;
				m_OOFValue = bestOOF;
			}
				
		}

		// new range
		std::transform (m_vRange.cbegin (), m_vRange.cend (), m_vRange.begin (), 
			[&](auto& pr) { return std::make_pair ((1.0-m_rScale)*pr.first, (1.0 - m_rScale)*pr.second); });

	} while (m_OOFValue > m_rTolerance);

	return;
}

//////////////////////////////////////////////////
