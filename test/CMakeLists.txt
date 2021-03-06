if(NOT BUILD_TESTING)
  return()
endif()

function(test_exe exename srcname)
  add_executable(${exename} ${srcname})
  target_link_libraries(${exename} PRIVATE GOAL)
endfunction()

function(mpi_test testname exe np)
  add_test(
    NAME ${testname}
    COMMAND ${MPIEXE} ${MPIFLAGS} ${np} "./${exe}" ${ARGN})
endfunction()

function(primal_test input)
  set(inyaml "${input}.yaml")
  set(exe ${CMAKE_CURRENT_BINARY_DIR}/../src/GoalPrimal)
  copy(${inyaml})
  add_test(${input} ${MPIEXE} ${MPIFLAGS} 4 ${exe} ${inyaml})
endfunction()

copy_dir(mesh)

set(square_1p_args
  "./mesh/square/square.dmg"
  "./mesh/square/square-serial.smb"
  "./mesh/square/square.txt")

set(square_4p_args
  "./mesh/square/square.dmg"
  "./mesh/square/square.smb"
  "./mesh/square/square.txt")

set(cube_1p_args
  "./mesh/cube/cube.dmg"
  "./mesh/cube/cube-serial.smb"
  "./mesh/cube/cube.txt")

set(cube_4p_args
  "./mesh/cube/cube.dmg"
  "./mesh/cube/cube.smb"
  "./mesh/cube/cube.txt")


test_exe(test_control control.cpp)
mpi_test(control_1p test_control 1)
mpi_test(control_4p test_control 4)

test_exe(test_disc disc.cpp)
mpi_test(disc_2D_1p test_disc 1 ${square_1p_args})
mpi_test(disc_2D_4p test_disc 4 ${square_4p_args})
mpi_test(disc_3D_1p test_disc 1 ${cube_1p_args})
mpi_test(disc_3D_4p test_disc 4 ${cube_4p_args})

test_exe(test_nested nested.cpp)
mpi_test(nested_2D_1p test_nested 1 ${square_1p_args})
mpi_test(nested_2D_4p test_nested 4 ${square_4p_args})
mpi_test(nested_3D_1p test_nested 1 ${cube_1p_args})
mpi_test(nested_3D_4p test_nested 4 ${cube_4p_args})

test_exe(test_sol_info sol_info.cpp)
mpi_test(sol_info_2D_1p test_sol_info 1 ${square_1p_args})
mpi_test(sol_info_2D_4p test_sol_info 4 ${square_4p_args})
mpi_test(sol_info_3D_1p test_sol_info 1 ${cube_1p_args})
mpi_test(sol_info_3D_4p test_sol_info 4 ${cube_4p_args})

test_exe(test_states states.cpp)
mpi_test(states_2D test_states 1 ${square_1p_args})
mpi_test(states_3D test_states 1 ${cube_1p_args})

test_exe(test_dofs dofs.cpp)
mpi_test(dofs_2D test_dofs 1 ${square_1p_args})
mpi_test(dofs_3D test_dofs 1 ${cube_1p_args})

bob_end_subdir()
