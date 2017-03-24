if (NOT BUILD_TESTING)
  return()
endif()

function(add_step number)
  set(exename step${number})
  set(srcname step${number}.cpp)
  add_executable(${exename} ${srcname})
  target_link_libraries(${exename} Goal)
  add_test(${exename} ${exename})
endfunction(add_step)

add_step(1)
add_step(2)
