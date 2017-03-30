if(NOT Goal_ENABLE_ELASTICITY)
  return()
endif()

message(STATUS "Setting up elasticity mini-app")

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/elasticity.dox
  ${CMAKE_CURRENT_BINARY_DIR}/elasticity.dox)

set(Goal_ELASTICITY_SOURCES
elast_physics.cpp
)

add_library(GoalElasticity ${Goal_ELASTICITY_SOURCES})
target_link_libraries(GoalElasticity PUBLIC Goal)

bob_end_subdir()
