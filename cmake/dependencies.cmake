# find SCOREC in the user-specified location
set(Goal_USE_SCOREC_DEFAULT ON)
bob_public_dep(SCOREC)

# find Trilinos in the user-specified location
set(Goal_USE_Trilinos_DEFAULT ON)
bob_public_dep(Trilinos)

# assert packages are enabled
macro(assert_trilinos_pkg pkg_name)
  list(FIND Trilinos_PACKAGE_LIST ${pkg_name} idx)
  if(NOT idx GREATER -1)
    message(FATAL_ERROR "Trilinos: ${pkg_name} not enabled")
  endif()
  message(STATUS "Trilinos: ${pkg_name} is enabled!")
endmacro()
assert_trilinos_pkg(Pamgen)
assert_trilinos_pkg(Phalanx)
assert_trilinos_pkg(Belos)
assert_trilinos_pkg(Ifpack2)
assert_trilinos_pkg(MiniTensor)
