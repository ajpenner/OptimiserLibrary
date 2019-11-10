// LineSearchFactory.h
//
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "stdafx.h"
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