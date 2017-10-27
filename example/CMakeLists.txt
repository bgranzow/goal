if(NOT BUILD_TESTING)
  return()
endif()

copy_dir(mesh)

add_subdirectory(primal)
add_subdirectory(spr)

bob_end_subdir()
