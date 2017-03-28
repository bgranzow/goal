if (NOT BUILD_TESTING)
  return()
endif()

set(MPIEXE mpirun)
set(MPIFLAGS -np)

function(copy INPUT)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/${INPUT}
    ${CMAKE_CURRENT_BINARY_DIR}/${INPUT})
endfunction()

function(add_serial_test EXE INPUT)
  set(INPUTYAML "in_${INPUT}.yaml")
  copy(${INPUTYAML})
  add_test(${INPUT} ${EXE} ${INPUTYAML})
endfunction()

function(add_parallel_test EXE INPUT)
  set(INPUTYAML "in_${INPUT}.yaml")
  copy(${INPUTYAML})
  add_test(${INPUT} ${MPIEXE} ${MPIFLAGS} 4 ${EXE} ${INPUTYAML})
endfunction()

add_subdirectory(steps)
add_subdirectory(poisson)

bob_end_subdir()
