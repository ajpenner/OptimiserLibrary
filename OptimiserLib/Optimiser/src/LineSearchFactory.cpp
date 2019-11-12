// LineSearchFactory.cpp
//
// 2017
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif
#include "LineSearchFactory.h"

//////////////////////////////////////////////////

std::unique_ptr<ILineSearch> LineSearchFactory::make_lineSearch (eLineSearch::type choice)
{
	if (choice == eLineSearch::eStrongWolfe)
		return std::make_unique<CLineSearchStrongWolfe> ();
	else if (choice == eLineSearch::eWeakWolfe)
		return std::make_unique<CLineSearchWeakWolfe> ();
	else if (choice == eLineSearch::eBackTracking)
		return std::make_unique<CLineSearchBackTrack> ();
	else if (choice == eLineSearch::eBackTracking)
		return std::make_unique<CLineSearchBackTrack> ();
	else
		return nullptr;
}

//////////////////////////////////////////////////
