if(NOT ${Goal_ENABLE_POISSON})
  return()
endif()

if("${Goal_DATA}" STREQUAL "")
  message(STATUS "WARNING: BUILD_TESTING=ON & Goal_ENABLE_POISSON=ON")
  message(STATUS "-> but Goal_DATA was not set")
  message(STATUS "-> skipping poisson tests")
  return()
endif()

message(STATUS "Setting up poisson mini-app tests")
set(POISSON_EXE ${CMAKE_BINARY_DIR}/miniapps/poisson/GoalPoissonSolver)
add_serial_test(${POISSON_EXE} poisson_lshape_serial)

bob_end_subdir()
