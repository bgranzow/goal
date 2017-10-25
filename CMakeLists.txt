cmake_minimum_required(VERSION 3.0.0)

project(GOAL VERSION 0.1.0 LANGUAGES CXX)

include(cmake/bob.cmake)
bob_begin_package()
bob_begin_cxx_flags()
bob_cxx11_flags()
bob_end_cxx_flags()

set(GOAL_FAD_SIZE "16" CACHE STRING "Maximum Sacado derivative array size")
option(GOAL_ENABLE_SNAPPING "Enable snapping via Simmetrix" OFF)
option(GOAL_ENABLE_MECH "Enable mechanics research code" OFF)
message(STATUS "GOAL_FAD_SIZE: ${GOAL_FAD_SIZE}")
message(STATUS "GOAL_ENABLE_SNAPPING: ${GOAL_ENABLE_SNAPPING}")

include(cmake/dependencies.cmake)
include(cmake/functions.cmake)

add_subdirectory(src)
add_subdirectory(test)

bob_end_package()
