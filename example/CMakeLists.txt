if(NOT BUILD_TESTING)
  return()
endif()

copy_dir(mesh)

add_subdirectory(primal)

bob_end_subdir()
