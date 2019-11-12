cmake_minimum_required(VERSION 3.0)

set(COMPONENT_VERSION_MAJOR 1)
set(COMPONENT_VERSION_MINOR 0)

set(COMPONENT_NAME "OptimiserFunctions")
set(COMPONENT_SOURCE_DIR ${PROJECT_SOURCE_DIR}/Optimiser/src)

#set(CMAKE_INSTALL_RPATH "/usr/lib/domino")
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_library(${COMPONENT_NAME} 
            SHARED
            ${COMPONENT_SOURCE_DIR}/Deb.cpp
            ${COMPONENT_SOURCE_DIR}/CompoundMOM.cpp
            ${COMPONENT_SOURCE_DIR}/CompoundPenalty.cpp
            ${COMPONENT_SOURCE_DIR}/Himmelblau.cpp
            ${COMPONENT_SOURCE_DIR}/LinearFunction.cpp
            ${COMPONENT_SOURCE_DIR}/Rosenbrock.cpp
            ${COMPONENT_SOURCE_DIR}/FirstConstraint.cpp
            ${COMPONENT_SOURCE_DIR}/SecondConstraint.cpp
            ${COMPONENT_SOURCE_DIR}/NegativeFirstConstraint.cpp
            ${COMPONENT_SOURCE_DIR}/Parabolic.cpp
            )

target_include_directories(${COMPONENT_NAME}
							PRIVATE
						${PROJECT_SOURCE_DIR}/Optimiser/include	
						  )


target_link_libraries(${COMPONENT_NAME}
                     )

set_target_properties(${COMPONENT_NAME} PROPERTIES 
                      VERSION ${COMPONENT_VERSION_MAJOR}.${COMPONENT_VERSION_MINOR} 
                      SOVERSION ${COMPONENT_VERSION_MAJOR} 
                     )

install (TARGETS ${COMPONENT_NAME} 
         LIBRARY DESTINATION lib/domino
)

install (DIRECTORY ${PROJECT_SOURCE_DIR}/include/
         DESTINATION include/InterComm
         FILES_MATCHING PATTERN "*.h"
)
