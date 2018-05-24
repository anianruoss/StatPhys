cmake_minimum_required (VERSION 3.1)

project (MDAtom)

set (CMAKE_CXX_STANDARD 17)

# set default build type to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

set (headers
     lib/MDParameters.h
     lib/ParameterIO.h
     lib/ParameterValidityChecker.h
     lib/CoordinatesAndVelocitiesInitializer.h
     lib/Timer.h
     lib/TrajectoryFileWriter.h
     lib/MDSimulation.h
     lib/MDRun.h
     lib/MDRunOutput.h
     lib/InteractionCalculator.h
     lib/RadialDistribution.h
     lib/InstantaneousRadialDistribution.h
     lib/AveragedRadialDistribution.h
     lib/BinaryIO.h
     lib/PeriodicBoundaryConditions.h
     lib/CenterOfMassCalculator.h
     lib/RandomNumberGenerator.h
)

set (sources
     lib/MDParameters.cpp
     lib/ParameterIO.cpp
     lib/ParameterValidityChecker.cpp
     lib/CoordinatesAndVelocitiesInitializer.cpp
     lib/Timer.cpp
     lib/TrajectoryFileWriter.cpp
     lib/MDSimulation.cpp
     lib/MDRun.cpp
     lib/MDRunOutput.cpp
     lib/InteractionCalculator.cpp
     lib/InstantaneousRadialDistribution.cpp
     lib/AveragedRadialDistribution.cpp
     lib/BinaryIO.cpp
     lib/PeriodicBoundaryConditions.cpp
     lib/CenterOfMassCalculator.cpp
     lib/RandomNumberGenerator.cpp
)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -march=native")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

# base
add_library(mdatomLib ${sources} ${headers})
set_target_properties(mdatomLib PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(mdatomLib PUBLIC lib)

add_executable(mdatom main.cpp)
target_link_libraries(mdatom mdatomLib)

# harmonic
add_library(mdatomLibHarmonic ${sources} ${headers})
set_target_properties(mdatomLibHarmonic PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(mdatomLibHarmonic PUBLIC lib)
target_compile_definitions(mdatomLibHarmonic PUBLIC HARMONIC)

add_executable(mdatomHarmonic main.cpp)
target_link_libraries(mdatomHarmonic mdatomLibHarmonic)

# shake
add_library(mdatomLibShake ${sources} ${headers})
set_target_properties(mdatomLibShake PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(mdatomLibShake PUBLIC lib)
target_compile_definitions(mdatomLibShake PUBLIC SHAKE)

add_executable(mdatomShake main.cpp)
target_link_libraries(mdatomShake mdatomLibShake)

