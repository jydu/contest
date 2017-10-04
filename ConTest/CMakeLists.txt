# CMake script for Bio++ Program Suite
# Author: Julien Dutheil
# Created: 17/09/2009

ADD_EXECUTABLE(contest ConTest.cpp)
TARGET_LINK_LIBRARIES(contest ${LIBS})
SET_TARGET_PROPERTIES(contest PROPERTIES LINK_SEARCH_END_STATIC ${BUILD_STATIC})


# Install libs
INSTALL(TARGETS contest DESTINATION bin)
