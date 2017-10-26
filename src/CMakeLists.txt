set(GOAL_SOURCES
goal_control.cpp
goal_disc.cpp
goal_nested.cpp
)

set(GOAL_HEADERS
goal_scalar_types.hpp
goal_data_types.hpp
goal_control.hpp
goal_disc.hpp
goal_nested.hpp
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

bob_end_subdir()
