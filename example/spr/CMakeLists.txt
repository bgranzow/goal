function(spr_test input)
  set(inyaml "${input}.yaml")
  set(exe ${CMAKE_BINARY_DIR}/src/GoalSPR)
  copy(${inyaml})
  add_test(${input} ${MPIEXE} ${MPIFLAGS} 4 ${exe} ${inyaml})
endfunction()

primal_test(neohookean_notch_2D)

bob_end_subdir()
