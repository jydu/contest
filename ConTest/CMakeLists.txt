# CMake script for ConTest
# Author: Julien Dutheil
# Created: 17/09/2009

ADD_EXECUTABLE(contest ConTest.cpp)
SET (contest-targets
	contest
  )

FOREACH (target ${contest-targets})
  # Link (static or shared)
  IF (BUILD_STATIC)
    TARGET_LINK_LIBRARIES (${target} ${BPP_LIBS_STATIC})
    SET_TARGET_PROPERTIES (${target} PROPERTIES LINK_SEARCH_END_STATIC TRUE)
  ELSE (BUILD_STATIC)
    TARGET_LINK_LIBRARIES (${target} ${BPP_LIBS_SHARED})
    SET_TARGET_PROPERTIES (${target} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
  ENDIF (BUILD_STATIC)
ENDFOREACH (target)

INSTALL (TARGETS ${contest-targets} DESTINATION ${CMAKE_INSTALL_BINDIR})
