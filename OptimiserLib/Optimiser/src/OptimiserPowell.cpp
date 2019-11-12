// OptimiserPowell.cpp : Defines the entry point for the console application.
//
// 2016
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif

//////////////////////////////////////////////////

#include "OptimiserPowell.h"
#include "MinimiserBisection.h" // need a line search in here
#include "MinimiserBrent.h" // need a line search in here
#include "MinimiserGolden.h" // need a line search in here
#include "optimisable.h"
#include "BracketBoundingPhase.h"
#include "LinearFunction.h"

//////////////////////////////////////////////////

void COptimiserPowell::SetParameters (const std::vector<variant>& parameters)
{
	assert (parameters.size () >= 1);
	m_alpha = boost::get<double>(parameters[ePowellParameters::eAlpha]);
}

//////////////////////////////////////////////////

double COptimiserPowell::GetMinimalScaleValue (std::shared_ptr<CLinearFunction> ptrVF, CMinimiserGolden& minimiser, CBracketBoundingPhase& range, const arma::vec& direction)
{
	ptrVF->SetVector (direction); // CLinearFunction needs to be renamed, perhaps CLinearFunction?
	ptrVF->SetOffset (m_MinVector);

	range.SetOptimisableFunction (ptrVF);
	range.SetInitialGuess (arma::vec({1.0}));
	range.Optimise(0); // this is broken
	auto values = range.GetRange();

	minimiser.SetOptimisableFunction (ptrVF);
	arma::vec guess = m_MinVector + direction;
	minimiser.SetInitialGuess (guess);
	minimiser.SetLeftBracket (arma::vec ({ values.at(0) }));
	minimiser.SetRightBracket (arma::vec ({ values.at(1) }));

	minimiser.Optimise (0);

	return minimiser.GetVector ().at (0);
}

//////////////////////////////////////////////////

void COptimiserPowell::Optimise ( double target )
{
	const auto iDim = m_ptrOptimisableFunction->GetFunctionDimension ();
	auto matEye = arma::mat (iDim, iDim).eye (); // fast way to obtain unit vectors

	// use a minimiser
	CMinimiserGolden minimiser;
	minimiser.SetOptimisableFunction (m_ptrOptimisableFunction);
	minimiser.SetTolerance(m_rTolerance);

	// calculate bounding phase
	CBracketBoundingPhase range;
	std::vector<variant> delta = { 1.0 };
	range.SetParameters (delta);

	auto ptrVF = std::make_shared<CLinearFunction>(m_ptrOptimisableFunction);

// below is messy
	std::vector<arma::vec> vVectorDirection(iDim);
	for (arma::uword i = 0; i < iDim; ++i) // replace with transform, really should just do a submat view
	{
		vVectorDirection[i] = matEye.col (i);
	}

	std::vector<arma::vec> coordinates(3);
	auto deltaOOF(std::numeric_limits<double>::max());
	do
	{
		for (arma::uword i = 0; i < iDim; ++i) // replace with transform
		{
			auto u = vVectorDirection[i];
			auto scale = GetMinimalScaleValue (ptrVF, minimiser, range, u);
			m_MinVector += scale*u;
			coordinates[i] = m_MinVector;
		}
		// along the first direction again
		auto u = vVectorDirection[0];
		auto scale = GetMinimalScaleValue (ptrVF, minimiser, range, u);
		m_MinVector += scale*u;
		deltaOOF = m_OOFValue;
		m_OOFValue = m_ptrOptimisableFunction->Evaluate(m_MinVector);
		deltaOOF -= m_OOFValue;
		coordinates[2] = m_MinVector;
		// swap direction vectors
		vVectorDirection[1] = vVectorDirection[0];
		vVectorDirection[0] = arma::normalise(coordinates[2] - coordinates[0]);

	} while (std::abs(m_OOFValue-target)>m_rTolerance && std::abs(deltaOOF) > m_rTolerance );

	return;
}

//////////////////////////////////////////////////
