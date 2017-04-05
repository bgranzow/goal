set(Goal_SOURCES
goal_control.cpp
goal_dimension.cpp
goal_discretization.cpp
goal_field.cpp
goal_indexer.cpp
goal_solution_info.cpp
goal_state_fields.cpp
goal_physics.cpp
goal_assembly.cpp
goal_output.cpp
goal_linear_solvers.cpp
goal_ev_gather_scalar.cpp
goal_ev_scalar_shape.cpp
goal_ev_interpolate_scalar.cpp
goal_ev_scatter_scalar.cpp
goal_ev_gather_vector.cpp
goal_ev_vector_shape.cpp
goal_ev_interpolate_vector.cpp
goal_ev_scatter_vector.cpp
goal_ev_scatter_functional.cpp
goal_ev_dirichlet_bcs.cpp
goal_ev_utils.cpp
)

set(Goal_HEADERS
goal_data_types.hpp
goal_control.hpp
goal_dimension.hpp
goal_workset.hpp
goal_discretization.hpp
goal_field_types.hpp
goal_field.hpp
goal_indexer.hpp
goal_solution_info.hpp
goal_state_fields.hpp
goal_physics.hpp
goal_traits.hpp
goal_assembly.hpp
goal_output.hpp
goal_linear_solvers.hpp
goal_ev_gather_scalar.hpp
goal_ev_scalar_shape.hpp
goal_ev_interpolate_scalar.hpp
goal_ev_scatter_scalar.hpp
goal_ev_gather_vector.hpp
goal_ev_vector_shape.hpp
goal_ev_interpolate_vector.hpp
goal_ev_scatter_vector.hpp
goal_ev_scatter_functional.hpp
goal_ev_dirichlet_bcs.hpp
goal_ev_utils.hpp
)

add_library(Goal ${Goal_SOURCES})
target_compile_definitions(Goal PUBLIC "-DGoal_FAD_SIZE=${Goal_FAD_SIZE}")
target_include_directories(Goal PUBLIC
  ${Trilinos_INCLUDE_DIRS}
  ${Trilinos_TPL_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(Goal PUBLIC
  ${Trilinos_LIBRARIES}
  ${Trilinos_TPL_LIBRARIES}
  ${Trilinos_EXTRA_LD_FLAGS})
if(SCOREC_FOUND)
  target_link_libraries(Goal PUBLIC SCOREC::core)
endif()
bob_export_target(Goal)

install(FILES ${Goal_HEADERS} DESTINATION include)

bob_end_subdir()
