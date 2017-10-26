function(primal_test input)
  set(inyaml "${input}.yaml")
  set(exe ${CMAKE_BINARY_DIR}/src/GoalPrimal)
  copy(${inyaml})
  add_test(${input} ${MPIEXE} ${MPIFLAGS} 4 ${exe} ${inyaml})
endfunction()

primal_test(neohookean_uniaxial_2D)
primal_test(neohookean_uniaxial_3D)
primal_test(J2_uniaxial_2D)
primal_test(J2_uniaxial_3D)

bob_end_subdir()
