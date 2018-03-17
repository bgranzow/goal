if(NOT BUILD_TESTING)
  return()
endif()

copy_dir(mesh)

add_subdirectory(primal)
add_subdirectory(restart)
add_subdirectory(spr)
add_subdirectory(nested)
add_subdirectory(manufactured)

bob_end_subdir()
