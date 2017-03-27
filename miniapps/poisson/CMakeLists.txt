if(NOT Goal_ENABLE_POISSON)
  return()
endif()

message(STATUS "Setting up poisson mini-app")

set(Goal_POISSON_SOURCES
poisson_ev_residual.cpp
)

add_library(GoalPoisson ${Goal_POISSON_SOURCES})
target_link_libraries(GoalPoisson PUBLIC Goal)
add_executable(GoalPoissonSolver poisson_solver.cpp)
target_link_libraries(GoalPoissonSolver GoalPoisson)
bob_export_target(GoalPoissonSolver)

bob_end_subdir()
