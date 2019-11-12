cmake_minimum_required(VERSION 3.0)

set(COMPONENT_VERSION_MAJOR 1)
set(COMPONENT_VERSION_MINOR 0)

set(COMPONENT_NAME "Optimisers")
set(COMPONENT_SOURCE_DIR ${PROJECT_SOURCE_DIR}/Optimiser/src)

#set(CMAKE_INSTALL_RPATH "/usr/lib/domino")
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_library(${COMPONENT_NAME} 
            SHARED
            ${COMPONENT_SOURCE_DIR}/BackTracking.cpp
            ${COMPONENT_SOURCE_DIR}/BracketBoundingPhase.cpp
            ${COMPONENT_SOURCE_DIR}/DirectionCauchy.cpp
            ${COMPONENT_SOURCE_DIR}/DirectionFactory.cpp
            ${COMPONENT_SOURCE_DIR}/DirectionNewton.cpp
            ${COMPONENT_SOURCE_DIR}/DirectionQuasiNewton.cpp
            ${COMPONENT_SOURCE_DIR}/ExtremumFactory.cpp
            ${COMPONENT_SOURCE_DIR}/LineSearchFactory.cpp
            ${COMPONENT_SOURCE_DIR}/MaximiserBrent.cpp
            ${COMPONENT_SOURCE_DIR}/MinimiserBisection.cpp
            ${COMPONENT_SOURCE_DIR}/MinimiserBoxed.cpp
            ${COMPONENT_SOURCE_DIR}/MinimiserBrent.cpp
            ${COMPONENT_SOURCE_DIR}/MinimiserExhaustive.cpp
            ${COMPONENT_SOURCE_DIR}/MinimiserGolden.cpp
            ${COMPONENT_SOURCE_DIR}/OptimiserBox.cpp
            ${COMPONENT_SOURCE_DIR}/OptimiserBrent.cpp
            ${COMPONENT_SOURCE_DIR}/OptimiserCauchy.cpp
            ${COMPONENT_SOURCE_DIR}/OptimiserCG.cpp
            ${COMPONENT_SOURCE_DIR}/OptimiserComplex.cpp
#            ${COMPONENT_SOURCE_DIR}/OptimiserEllipsoid.cpp
            ${COMPONENT_SOURCE_DIR}/OptimiserFactory.cpp
            ${COMPONENT_SOURCE_DIR}/OptimiserFactoryVM.cpp
#            ${COMPONENT_SOURCE_DIR}/OptimiserGRG.cpp
#            ${COMPONENT_SOURCE_DIR}/OptimiserHennig.cpp
#            ${COMPONENT_SOURCE_DIR}/OptimiserHJ.cpp
            ${COMPONENT_SOURCE_DIR}/OptimiserMOM.cpp
            ${COMPONENT_SOURCE_DIR}/OptimiserPenalty.cpp
            ${COMPONENT_SOURCE_DIR}/OptimiserRandom.cpp
            ${COMPONENT_SOURCE_DIR}/OptimiserVMLBFGS.cpp
            ${COMPONENT_SOURCE_DIR}/OptimiserMarquardt.cpp
            ${COMPONENT_SOURCE_DIR}/OptimiserNewton.cpp
            ${COMPONENT_SOURCE_DIR}/OptimiserNonParHen.cpp
            ${COMPONENT_SOURCE_DIR}/OptimiserPowell.cpp
            ${COMPONENT_SOURCE_DIR}/OptimiserSimplex.cpp
            ${COMPONENT_SOURCE_DIR}/OptimiserVMBFGS.cpp
            ${COMPONENT_SOURCE_DIR}/OptimiserVMBlend.cpp
            ${COMPONENT_SOURCE_DIR}/OptimiserVMBroyden.cpp
            ${COMPONENT_SOURCE_DIR}/OptimiserVMDFP.cpp
            ${COMPONENT_SOURCE_DIR}/OptimiserVMPSB.cpp
            ${COMPONENT_SOURCE_DIR}/OptimiserVMSR1.cpp
            ${COMPONENT_SOURCE_DIR}/StrongWolfe.cpp
            ${COMPONENT_SOURCE_DIR}/WeakWolfe.cpp
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
