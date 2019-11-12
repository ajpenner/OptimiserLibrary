// DirectionFactory.cpp
//
// 2017
//////////////////////////////////////////////////

#ifdef _WIN32
#include "stdafx.h"
#endif
#include "DirectionFactory.h"

//////////////////////////////////////////////////

std::unique_ptr<IDirection> DirectionFactory::make_direction (eDirection::type choice)
{
	if (choice == eDirection::eCauchy)
		return std::make_unique<CDirectionCauchy> ();
	else if (choice == eDirection::eNewton)
		return std::make_unique<CDirectionNewton> ();
	else if (choice == eDirection::eQuasiNewton)
		return std::make_unique<CDirectionQuasiNewton> ();
	else
		return nullptr;
}

//////////////////////////////////////////////////
