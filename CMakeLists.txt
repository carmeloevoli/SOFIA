project(ifeda CXX)

cmake_minimum_required(VERSION 2.6)

set( CMAKE_VERBOSE_MAKEFILE off )
set( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" )
set( CMAKE_CXX_COMPILER "g++-9" )
set( CMAKE_CXX_FLAGS "-std=c++11 -O3" )
set( IFEDA_EXTRA_INCLUDES )
set( IFEDA_EXTRA_LIBRARIES )

include_directories( "${CMAKE_SOURCE_DIR}" )

set(IFEDA_SRCS 
	src/advection.cpp
	src/axis.cpp
	src/diffusion.cpp
	src/evolutor.cpp
	src/losses.cpp
	src/main.cpp
	src/output.cpp
	src/params.cpp
	src/particle.cpp
	src/sources.cpp
	src/tridiag.cpp
	src/utils.cpp
	src/waves.cpp
    )

# GSL (external)
find_package(GSL)
if(GSL_FOUND)
	set(HAVE_LIBGSL 1)
	list(APPEND IFEDA_EXTRA_INCLUDES ${GSL_INCLUDE_DIR})
	list(APPEND IFEDA_EXTRA_LIBRARIES ${GSL_LIBRARIES})
else(GSL_FOUND)
	message(FATAL_ERROR "Require GSL. Set GSL_HOME")
endif(GSL_FOUND)

# MAKE EXECUTABLE
include_directories (include ${IFEDA_EXTRA_INCLUDES})

add_executable(ifeda ${IFEDA_SRCS} )
target_link_libraries (ifeda ${IFEDA_EXTRA_LIBRARIES})

# Show summary.
message (STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
