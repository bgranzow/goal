set(GOAL_SOURCES
goal_control.cpp
goal_disc.cpp
goal_nested.cpp
goal_sol_info.cpp
goal_integrator.cpp
goal_soln.cpp
goal_scalar_weight.cpp
goal_soln_adjoint.cpp
goal_residual.cpp
goal_qoi.cpp
goal_point_wise.cpp
goal_assembly.cpp
goal_dbcs.cpp
goal_linear_solve.cpp
goal_physics.cpp
goal_primal.cpp
goal_functional.cpp
goal_nested_adjoint.cpp
goal_error.cpp
goal_size_field.cpp
goal_output.cpp
)

set(GOAL_HEADERS
goal_scalar_types.hpp
goal_data_types.hpp
goal_eval_modes.hpp
goal_control.hpp
goal_disc.hpp
goal_nested.hpp
goal_sol_info.hpp
goal_integrator.hpp
goal_soln.hpp
goal_scalar_weight.hpp
goal_soln_adjoint.hpp
goal_residual.hpp
goal_qoi.hpp
goal_point_wise.hpp
goal_assembly.hpp
goal_dbcs.hpp
goal_linear_solve.hpp
goal_physics.hpp
goal_primal.hpp
goal_functional.hpp
goal_nested_adjoint.hpp
goal_error.hpp
goal_size_field.hpp
goal_output.hpp
)

add_library(GOAL ${GOAL_SOURCES} ${GOAL_HEADERS})
target_compile_definitions(GOAL PUBLIC
  "-DGOAL_FAD_SIZE=${GOAL_FAD_SIZE}")
if(GOAL_ENABLE_SNAPPING)
  target_compile_definitions(GOAL PUBLIC
    "-DGOAL_ENABLE_SNAPPING=1")
endif()
target_include_directories(GOAL PUBLIC
  ${Trilinos_INCLUDE_DIRS}
  ${Trilinos_TPL_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(GOAL PUBLIC
  ${Trilinos_LIBRARIES}
  ${Trilinos_TPL_LIBRARIES}
  ${Trilinos_EXTRA_LD_FLAGS}
  SCOREC::core)
bob_export_target(GOAL)

function(add_exe exename exesrc)
  add_executable(${exename} ${exesrc})
  target_link_libraries(${exename} PRIVATE GOAL)
  bob_export_target(${exename})
endfunction()

add_exe(GoalPrimal main_primal.cpp)
add_exe(GoalAdjoint main_nested.cpp)

bob_end_subdir()
