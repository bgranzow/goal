cmake_minimum_required(VERSION 3.0.0)

project(Goal VERSION 0.1.0 LANGUAGES CXX)

include(cmake/bob.cmake)

bob_begin_package()

set(Scorec_PREFIX "" CACHE PATH "Scorec installation directory")
set(Trilinos_PREFIX "" CACHE PATH "Trilinos installation directory")
set(GOAL_FAD_SIZE "40" CACHE STRING "Maximum Sacado derivative array size")
option(GOAL_ENABLE_POISSON "Enable the poisson's equation mini-app" OFF)
option(GOAL_ENABLE_ELASTICITY "Enable the elasticity mini-app" OFF)
option(GOAL_ENABLE_ALL_APPS "Enable all of the mini-apps provided by Goal" OFF)
if(GOAL_ENABLE_ALL_APPS)
  set(GOAL_ENABLE_POISSON ON)
  set(GOAL_ENABLE_ELASTICITY ON)
endif()

bob_begin_cxx_flags()
bob_cxx11_flags()
bob_end_cxx_flags()

message(STATUS "GOAL_FAD_SIZE: ${GOAL_FAD_SIZE}")
message(STATUS "GOAL_ENABLE_POISSON: ${GOAL_ENABLE_POISSON}")
message(STATUS "GOAL_ENABLE_ELASTICITY: ${GOAL_ENABLE_ELASTICITY}")
message(STATUS "GOAL_ENABLE_ALL_APPS: ${GOAL_ENABLE_ALL_APPS}")

include(cmake/dependencies.cmake)

add_subdirectory(src)
add_subdirectory(miniapps)
add_subdirectory(test)
add_subdirectory(doc)

bob_end_package()
