function(nested_test input)
  set(inyaml "${input}.yaml")
  set(exe ${CMAKE_BINARY_DIR}/src/GoalNested)
  copy(${inyaml})
  add_test(${input} ${MPIEXE} ${MPIFLAGS} 4 ${exe} ${inyaml})
endfunction()

nested_test(cooks_point_wise)
nested_test(cooks_avg_disp_subdomain)
nested_test(notch_2D_max_vm)

bob_end_subdir()
