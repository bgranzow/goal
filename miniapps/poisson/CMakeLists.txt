if(NOT Goal_ENABLE_POISSON)
  return()
endif()

message(STATUS "Setting up poisson mini-app")

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/poisson.dox
  ${CMAKE_CURRENT_BINARY_DIR}/poisson.dox)

set(Goal_POISSON_SOURCES
poisson_ev_residual.cpp
poisson_physics.cpp
)

add_library(GoalPoisson ${Goal_POISSON_SOURCES})
target_link_libraries(GoalPoisson PUBLIC Goal)
add_executable(GoalPoissonSolver poisson_solver.cpp)
target_link_libraries(GoalPoissonSolver PRIVATE GoalPoisson)
bob_export_target(GoalPoissonSolver)

bob_end_subdir()
