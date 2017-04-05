find_package(Doxygen)
if(NOT DOXYGEN_FOUND)
  return()
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
  ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/meshgen.dox
  ${CMAKE_CURRENT_BINARY_DIR}/meshgen.dox)

add_custom_target(doc
  ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating API documentation with Doxygen" VERBATIM)
