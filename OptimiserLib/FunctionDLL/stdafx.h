// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
// Windows Header Files:

#define NOMINMAX // remove windows default min & max definitions
#include <windows.h>

// global definitions
#define realEmpty std::numeric_limits<double>::max()
#define realZero 0x000000000

// general includes
#include <assert.h>