cmake_minimum_required(VERSION 3.0.0)

project(Goal VERSION 0.1.0 LANGUAGES CXX)

include(cmake/bob.cmake)

bob_begin_package()

set(Scorec_PREFIX "" CACHE PATH "Scorec installation directory")
set(Trilinos_PREFIX "" CACHE PATH "Trilinos installation directory")
set(GOAL_FAD_SIZE "40" CACHE STRING "Maximum Sacado derivative array size")

bob_begin_cxx_flags()
bob_cxx11_flags()
bob_end_cxx_flags()

message(STATUS "GOAL_FAD_SIZE: ${GOAL_FAD_SIZE}")

include(cmake/dependencies.cmake)

add_subdirectory(src)
add_subdirectory(doc)
add_subdirectory(test)

bob_end_package()
