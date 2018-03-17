function(mfct_primal_test input)
  set(inyaml "${input}.yaml")
  set(exe ${CMAKE_BINARY_DIR}/src/GoalPrimal)
  copy(${inyaml})
  add_test(${input} ${MPIEXE} ${MPIFLAGS} 1 ${exe} ${inyaml})
endfunction()

mfct_primal_test(elastic_primal_2D)

bob_end_subdir()
