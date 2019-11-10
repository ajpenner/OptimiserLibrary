// ExtremumFactory.h
//
// 2017
//////////////////////////////////////////////////

#pragma once

//////////////////////////////////////////////////

#include "stdafx.h"
#include "ExtremumFactory.h"

//////////////////////////////////////////////////

std::unique_ptr<BExtremum> ExtremumFactory::make_extremum (eExtremum::type choice)
{
	if (choice == eExtremum::eMaxBrent)
		return std::make_unique<CMaximiserBrent> ();
	else if (choice == eExtremum::eMinBrent)
		return std::make_unique<CMinimiserBrent> ();
	else if (choice == eExtremum::eBisection)
		return std::make_unique<CMinimiserBisection> ();
	else if (choice == eExtremum::eBoxed)
		return std::make_unique<CMinimiserBoxed> ();
	else if (choice == eExtremum::eExhaustive)
		return std::make_unique<CMinimiserExhaustive> ();
	else if (choice == eExtremum::eGolden)
		return std::make_unique<CMinimiserGolden> ();
	else
		return nullptr;
}

//////////////////////////////////////////////////