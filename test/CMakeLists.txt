if(NOT BUILD_TESTING)
  return()
endif()

set(MPIEXE mpirun)
set(MPIFLAGS -np)

function(test_exe exename srcname)
  add_executable(${exename} ${srcname})
  target_link_libraries(${exename} PRIVATE Goal)
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

test_exe(disc disc.cpp)
mpi_test(disc2D_1p disc 1 "box2D.dmg" "box2D_1p.smb" "box2D.txt")
mpi_test(disc2D_4p disc 4 "box2D.dmg" "box2D_4p.smb" "box2D.txt")

add_custom_target(pretest COMMAND)
add_dependencies(pretest meshgen)
