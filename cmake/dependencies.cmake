# find trilinos only in the user-specified location
find_package(Trilinos 12.0.0 REQUIRED PATHS ${Trilinos_PREFIX} NO_DEFAULT_PATH)

# purge duplicates in trilinos lists
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

# ensure needed trilinos packages were built
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

# find SCOREC in the user-specified location
set(Goal_USE_SCOREC_DEFAULT ON)
bob_public_dep(SCOREC)

# print the locations of the found packages
message(STATUS "found Trilinos: ${Trilinos_DIR} (${Trilinos_VERSION})")
message(STATUS "found SCOREC: ${SCOREC_DIR} (${SCOREC_VERSION})")
