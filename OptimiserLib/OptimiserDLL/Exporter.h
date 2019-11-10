//
// Define an export method for use with the optimiser DLL
//

#pragma once

#if !defined(OPTIMISERDLL_EXPORTS)
#   pragma comment(lib, "OptimiserDLL.lib")
#endif

#ifdef OPTIMISERDLL_EXPORTS
#define OPTIMISER __declspec(dllexport)
#else
#define OPTIMISER __declspec(dllimport)
#endif

#if !defined(FUNCTIONDLL_EXPORTS)
#   pragma comment(lib, "FunctionDLL.lib")
#endif

#ifdef FUNCTIONDLL_EXPORTS
#define FUNCTION __declspec(dllexport)
#else
#define FUNCTION __declspec(dllimport)
#endif