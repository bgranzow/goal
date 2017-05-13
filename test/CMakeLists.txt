if(NOT BUILD_TESTING)
  return()
endif()

set(MPIEXE mpirun)
set(MPIFLAGS -np)

function(test_exe exename srcname)
  add_executable(${exename} ${srcname})
  target_link_libraries(${exename} PRIVATE Goal)
  add_dependencies(${exename} box2D_1p)
  add_dependencies(${exename} box2D_4p)
  add_dependencies(${exename} box3D_1p)
  add_dependencies(${exename} box3D_4p)
  add_dependencies(${exename} boxAssoc)
endfunction()

function(mpi_test testname exe np)
  add_test(
    NAME ${testname}
    COMMAND ${MPIEXE} ${MPIFLAGS} ${np} "./${exe}" ${ARGN})
endfunction()

include(meshgen.cmake)

test_exe(control control.cpp)
mpi_test(control_1p control 1)
mpi_test(control_4p control 4)
