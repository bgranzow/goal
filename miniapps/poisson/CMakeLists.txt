if(NOT GOAL_ENABLE_POISSON)
  return()
endif()

message(STATUS "Setting up poisson mini-app")

set(GOAL_POISSON_SOURCES
poisson_ev_residual.cpp
)

add_library(GoalPoisson ${GOAL_POISSON_SOURCES})
target_link_libraries(GoalPoisson PUBLIC Goal)

bob_end_subdir()
