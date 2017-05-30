set(Goal_USE_SCOREC_DEFAULT ON)
set(Goal_USE_Trilinos_DEFAULT ON)
bob_public_dep(Trilinos)
bob_public_dep(SCOREC)

function(assert_trilinos_pkg pkg_name)
  list(FIND Trilinos_PACKAGE_LIST ${pkg_name} idx)
  if(NOT idx GREATER -1)
    message(FATAL_ERROR "Trilinos: ${pkg_name} not enabled")
  endif()
  message(STATUS "Trilinos: ${pkg_name} is enabled!")
endfunction()

assert_trilinos_pkg(Pamgen)
assert_trilinos_pkg(Phalanx)
assert_trilinos_pkg(Belos)
assert_trilinos_pkg(Ifpack2)
assert_trilinos_pkg(MiniTensor)
assert_trilinos_pkg(Amesos2)
assert_trilinos_pkg(MueLu)

list(REVERSE Trilinos_INCLUDE_DIRS)
list(REVERSE Trilinos_TPL_INCLUDE_DIRS)
list(REVERSE Trilinos_LIBRARIES)
list(REVERSE Trilinos_TPL_LIBRARIES)
list(REMOVE_DUPLICATES Trilinos_INCLUDE_DIRS)
list(REMOVE_DUPLICATES Trilinos_TPL_INCLUDE_DIRS)
list(REMOVE_DUPLICATES Trilinos_LIBRARIES)
list(REMOVE_DUPLICATES Trilinos_TPL_LIBRARIES)
list(REVERSE Trilinos_INCLUDE_DIRS)
list(REVERSE Trilinos_TPL_INCLUDE_DIRS)
list(REVERSE Trilinos_LIBRARIES)
list(REVERSE Trilinos_TPL_LIBRARIES)
