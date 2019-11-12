cmake_minimum_required(VERSION 3.0)

project(optimiser)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CONFIGURATION_TYPES "")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -s")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -s")
set(CMAKE_CXX_FLAGS_DEBUG  "-g")

add_definitions("-DBLOCKING")

set(CMAKE_DEBUG_POSTFIX "-dbg")

find_package(Boost REQUIRED)
find_package(Armadillo REQUIRED)

add_subdirectory(FunctionDLL)
add_subdirectory(OptimiserDLL)
#add_subdirectory(OptimiserTestDrive)

