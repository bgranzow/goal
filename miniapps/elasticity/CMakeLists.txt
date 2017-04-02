if(NOT Goal_ENABLE_ELASTICITY)
  return()
endif()

message(STATUS "Setting up elasticity mini-app")

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/elasticity.dox
  ${CMAKE_CURRENT_BINARY_DIR}/elasticity.dox)

set(Goal_ELASTICITY_SOURCES
elast_ev_stress.cpp
elast_ev_residual.cpp
elast_physics.cpp
)

add_library(GoalElasticity ${Goal_ELASTICITY_SOURCES})
target_link_libraries(GoalElasticity PUBLIC Goal)
add_executable(GoalElasticitySolver elast_solver.cpp)
target_link_libraries(GoalElasticitySolver PRIVATE GoalElasticity)
bob_export_target(GoalElasticitySolver)

bob_end_subdir()
