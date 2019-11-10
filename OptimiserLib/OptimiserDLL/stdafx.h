// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>

// STL
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <numeric>

// boost
#include <boost/iterator/zip_iterator.hpp>

// global definitions
#define realEmpty std::numeric_limits<double>::max()
#define realZero 0x000000000
