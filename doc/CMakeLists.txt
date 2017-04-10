find_package(Doxygen)
if(NOT DOXYGEN_FOUND)
  return()
endif()

if("${Goal_DATA}" STREQUAL "")
  message(STATUS "skipping documentation...")
  message(STATUS "-> because Goal_DATA was not set")
  return()
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
  ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

function(doc_config fname)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/${fname}
    ${CMAKE_CURRENT_BINARY_DIR}/${fname} @ONLY)
endfunction()

doc_config(meshgen.dox)

add_custom_target(doc
  ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating API documentation with Doxygen" VERBATIM)
