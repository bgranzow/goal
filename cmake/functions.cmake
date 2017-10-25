set(MPIEXE mpirun)
set(MPIFLAGS -np)

function(copy input)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/${input}
    ${CMAKE_CURRENT_BINARY_DIR}/${input} COPYONLY)
endfunction()

function(copy_dir input)
  file(
    COPY ${CMAKE_CURRENT_SOURCE_DIR}/${input}
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endfunction()
