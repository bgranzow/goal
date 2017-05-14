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

set(2D1p_args "box2D.dmg" "box2D_1p.smb" "box2D.txt")
set(2D4p_args "box2D.dmg" "box2D_4p.smb" "box2D.txt")
set(3D1p_args "box3D.dmg" "box3D_1p.smb" "box3D.txt")
set(3D4p_args "box3D.dmg" "box3D_4p.smb" "box3D.txt")

test_exe(control control.cpp)
mpi_test(control_1p control 1)
mpi_test(control_4p control 4)

test_exe(disc disc.cpp)
mpi_test(disc2D_1p disc 1 ${2D1p_args})
mpi_test(disc2D_4p disc 4 ${2D4p_args})
mpi_test(disc3D_4p disc 4 ${3D4p_args})

test_exe(field field.cpp)
mpi_test(field2D field 1 ${2D1p_args})
mpi_test(field3D field 1 ${3D1p_args})

add_custom_target(pretest COMMAND)
add_dependencies(pretest meshgen)
