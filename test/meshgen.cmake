function(find_scorec_exe exename)
  find_program(${exename}_exe ${exename}
    PATHS "${SCOREC_PREFIX}/bin" NO_DEFAULT_PATH)
  if(NOT ${exename_exe})
    message(FATAL_ERROR "scorec utility: ${exename} not found")
  endif()
endfunction()

find_scorec_exe(box)
find_scorec_exe(split)
message(STATUS "box exe: ${box_exe}")
message(STATUS "split exe: ${split_exe}")

add_custom_target(box2D_1p
  COMMAND ${box_exe}
  "10" "10" "0" "1" "1" "0" "1" "box2D.dmg" "box2D_1p.smb"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating serial 2D box" VERBATIM)

add_custom_target(box3D_1p
  COMMAND ${box_exe}
  "5" "5" "5" "1" "1" "1" "1" "box3D.dmg" "box3D_1p.smb"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating serial 3D box" VERBATIM)

add_custom_target(box2D_4p
  COMMAND ${MPIEXE} ${MPIFLAGS} 4 ${split_exe}
  "box2D.dmg" "box2D_1p.smb" "box2D_4p.smb" "4"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating parallel 2D box" VERBATIM)
add_dependencies(box2D_4p box2D_1p)

add_custom_target(box3D_4p
  COMMAND ${MPIEXE} ${MPIFLAGS} 4 ${split_exe}
  "box3D.dmg" "box3D_1p.smb" "box3D_4p.smb" "4"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating parallel 3D box" VERBATIM)
add_dependencies(box3D_4p box3D_1p)

add_executable(box_assoc box_assoc.cpp)
add_custom_target(boxAssoc
  COMMAND box_assoc
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating box associations" VERBATIM)

add_custom_target(meshgen
  COMMAND make
  "box2D_1p" "box2D_4p" "box3D_1p" "box3D_4p" "boxAssoc")
