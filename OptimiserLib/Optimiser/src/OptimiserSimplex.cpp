// COptimiserSimplex.cpp
//
// 2016
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif

//////////////////////////////////////////////////

#include "OptimiserSimplex.h"
#include "optimisable.h"
#include <cassert>

//////////////////////////////////////////////////

COptimiserSimplex::COptimiserSimplex () :
	m_gamma(1.5),
	m_beta(0.5),
	m_delta(0)
{ 
}

//////////////////////////////////////////////////

void COptimiserSimplex::SetLeftBracket (const arma::vec& vLeft)
{
	m_LeftBracket = vLeft; 
}

//////////////////////////////////////////////////

void COptimiserSimplex::SetRightBracket (const arma::vec& vRight)
{ 
	m_RightBracket = vRight; 
}

//////////////////////////////////////////////////

void COptimiserSimplex::SetParameters (const std::vector<variant>& parameters)
{
	assert (parameters.size () >= 2);
	m_beta = boost::get<double> (parameters[eSimplexParameters::eBeta]);
	m_gamma = boost::get<double> (parameters[eSimplexParameters::eGamma]);

	m_simplex.clear ();
	for (auto it = parameters.cbegin () + eSimplexParameters::eSimplexStart;
		it != parameters.cend (); ++it)
	{
		const auto& v = boost::get<arma::vec>(*it);
		m_simplex.push_back(v);
	}
}

//////////////////////////////////////////////////

double COptimiserSimplex::GetInverseOOF() const
{
	return m_InverseOOFValue;
}

//////////////////////////////////////////////////

void COptimiserSimplex::InverseOptimise(double target)
{
	Optimise(-target);
}

//////////////////////////////////////////////////

void COptimiserSimplex::CreateSimplex() // I do not understand this yet
{
	// use the random start created in OptimiserComplex
}

//////////////////////////////////////////////////

void COptimiserSimplex::Optimise (double target)
{
	assert (m_ptrOptimisableFunction);
	if (m_ptrOptimisableFunction->GetFunctionDimension () != 3)
	{
		m_delta = (std::sqrt (m_ptrOptimisableFunction->GetFunctionDimension () + 1) - 2) / (m_ptrOptimisableFunction->GetFunctionDimension () - 3);
	}
	else
	{
		m_delta = 0.25;
	}

	std::vector<double> eval;
	for(const arma::vec& vec : m_simplex)
	{
		eval.push_back(m_ptrOptimisableFunction->Evaluate(vec));
	}

	auto terminate = std::numeric_limits<double>::max();
	do
	{
		auto itWorst = std::max_element (eval.cbegin (), eval.cend ());
		arma::vec xW = m_simplex[std::distance (eval.cbegin (), itWorst)];
		auto fW = *itWorst;
		m_simplex.erase (m_simplex.begin () + std::distance (eval.cbegin (), itWorst));
		eval.erase( std::remove(eval.begin(), eval.end(), *itWorst), eval.end ());

		auto itWorstRemaining = std::max_element(eval.begin(), eval.end());
		auto thing = std::distance (eval.begin (), itWorstRemaining);
		auto& xG = m_simplex[std::distance(eval.begin(), itWorstRemaining)];
		auto fG = *itWorstRemaining;

		auto itBest = std::min_element (eval.cbegin (), eval.cend ());
		auto& xB = m_simplex[std::distance (eval.cbegin (), itBest)];
		auto fB = *itBest;

		auto xC = arma::vec(m_ptrOptimisableFunction->GetFunctionDimension ());
		xC.zeros ();
		for (auto it = m_simplex.cbegin (); it != m_simplex.cend (); ++it)
		{
			xC += *it;
		}
		xC /= static_cast<double>(m_ptrOptimisableFunction->GetFunctionDimension ()); // dimension is a uint8
		auto fC = m_ptrOptimisableFunction->Evaluate (xC);

		auto xR = (2 * xC - xW).eval ();

		auto xNew = xR;

		auto fR = m_ptrOptimisableFunction->Evaluate (xR);
		if (fR < fB)
		{
			xNew = ((1 + m_gamma)*xC - m_gamma*xW).eval ();
		}
		else if (fR >= fW)
		{
			xNew = ((1 - m_beta)*xC + m_beta*xW).eval ();
		}
		else if (fG < fR && fR < fW)
		{
			xNew = ((1 + m_beta)*xC - m_beta*xW).eval ();
		}

		m_simplex.push_back( xNew );
		eval.push_back( m_ptrOptimisableFunction->Evaluate(xNew) );

		terminate = 0x0; // should use std::accumulate, but there is a bug using armadillo in a lambda
		for (auto it = eval.cbegin (); it != eval.cend (); ++it)
		{
			terminate += std::pow(*it-fC,2);
		}
		terminate /= (m_ptrOptimisableFunction->GetFunctionDimension () + 1);
		terminate = std::sqrt (terminate);
	} while (terminate > m_rTolerance);

	auto itBest = std::min_element (eval.cbegin (), eval.cend ());
	m_MinVector = m_simplex[std::distance (eval.cbegin (), itBest)];
	m_OOFValue = *itBest;

	return;
}

//////////////////////////////`////////////////////
