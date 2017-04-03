if(NOT ${Goal_ENABLE_ELASTICITY})
  return()
endif()

if("${Goal_DATA}" STREQUAL "")
  message(STATUS "WARNING: BUILD_TESTING=ON & Goal_ENABLE_ELASTICITY=ON")
  message(STATUS "-> but Goal_DATA was not set")
  message(STATUS "-> skipping elasticity tests")
  return()
endif()

message(STATUS "Setting up elasticity mini-app tests")
set(ELAST_EXE ${CMAKE_BINARY_DIR}/miniapps/elasticity/GoalElasticitySolver)
add_serial_test(${ELAST_EXE} elasticity_square_serial)
add_parallel_test(${ELAST_EXE} elasticity_square)

bob_end_subdir()
