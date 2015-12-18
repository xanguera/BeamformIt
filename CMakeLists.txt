################################################################################
# BeamformIt: acoustic synchronization and beamforming toolkit
################################################################################

cmake_minimum_required (VERSION 2.6)
project (beamformIt)

##where to put the executables
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

##general options
#set (CMAKE_VERBOSE_MAKEFILE    ON)
#set (CMAKE_BUILD_TYPE          DEBUG)
set (CMAKE_BUILD_TYPE          RELEASE)

##C and C++ flags
set (CMAKE_C_FLAGS_DEBUG           "-D_REENTRANT -Wall -DDEPURACION -O0 -g3 -Wno-unused-variable")
set (CMAKE_CXX_FLAGS_DEBUG         "-ggdb -DLINUX -Wall -W -Wno-unknown-pragmas -Wno-unused-variable")
set (CMAKE_C_FLAGS_RELEASE         "-D_REENTRANT -Wall -DDEPURACION -DNDEBUG -O3 -Wno-unused-variable")
set (CMAKE_CXX_FLAGS_RELEASE       "-DLINUX -Wall -W -O3 -DNDEBUG -Wno-unknown-pragmas -Wno-unused-result -Wno-unused-variable")

set (LIBSND_INSTALL_DIR "" CACHE FILEPATH "Directory where libsndfile was installed")

##find libraries installed in the system
FIND_LIBRARY(LIBSNDFILE sndfile ${LIBSND_INSTALL_DIR}/lib)

##set particular settings for each architecture
set (FFTREAL_HDRS ${CMAKE_CURRENT_SOURCE_DIR}/external/FFTReal-2.11)
if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  set(os_link_libs /usr/local/lib /opt/local/lib)
  set(os_include_libs /opt/local/include ${FFTREAL_HDRS})
  set(os_libs ${LIBSNDFILE})
elseif( ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
  set(os_link_libs /usr/lib)
  set(os_include_libs /opt/local/include ${FFTREAL_HDRS})
  set(os_libs ${LIBSNDFILE})
endif()

##set some generic link directories
LINK_DIRECTORIES(${os_link_libs})
set (CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem") #needed for the SYSTEM link directories to work on osx
INCLUDE_DIRECTORIES(SYSTEM ${os_include_libs} ${LIBSND_INSTALL_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/includes)

#gather project files (note that when adding a new file cmake needs to be rerun)
file(GLOB src_files src/*.cc)
file(GLOB inc_files includes/*.h)
#take out those files that contain the main() function
file(GLOB to_remove src/BeamformIt.cc)
list(REMOVE_ITEM src_files ${to_remove})

##executable definition
ADD_EXECUTABLE(BeamformIt src/BeamformIt ${src_files} ${inc_files})
TARGET_LINK_LIBRARIES(BeamformIt ${os_libs} )
SET_TARGET_PROPERTIES(BeamformIt PROPERTIES LINKER_LANGUAGE CXX)

##doxygen target to generate API documentation with Doxygen (run with "make documentation")
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cfg-files/doxygen.cfg ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(documentation
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
  )
endif(DOXYGEN_FOUND)
