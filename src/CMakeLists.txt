set(Goal_SOURCES
goal_control.cpp
goal_dimension.cpp
goal_discretization.cpp
goal_field.cpp
goal_indexer.cpp
goal_solution_info.cpp
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
goal_traits.hpp
)

add_library(Goal ${Goal_SOURCES})
target_compile_definitions(Goal PUBLIC "-DGOAL_FAD_SIZE=${GOAL_FAD_SIZE}")
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
