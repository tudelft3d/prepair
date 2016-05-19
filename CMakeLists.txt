# Copyright (c) 2009-2014,
# Ken Arroyo Ohori    g.a.k.arroyoohori@tudelft.nl
# Hugo Ledoux         h.ledoux@tudelft.nl
# Martijn Meijers     b.m.meijers@tudelft.nl
# All rights reserved.
# 
# This file is part of prepair: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Licensees holding a valid commercial license may use this file in
# accordance with the commercial license agreement provided with
# the software.
# 
# This file is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#

project( prepair )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

# CGAL
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )
  message(SEND_ERROR "prepair requires the CGAL library")
  return()  
endif()

# include helper file
include( ${CGAL_USE_FILE} )

# Boost
find_package( Boost REQUIRED COMPONENTS program_options)

if ( NOT Boost_FOUND )
  message(SEND_ERROR "prepair requires the Boost library")
  return()  
endif()

INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

# GDAL
find_package( GDAL )

if ( NOT GDAL_FOUND )
  message(SEND_ERROR "prepair requires the GDAL library")
endif()

include_directories( ${GDAL_INCLUDE_DIR} )

# Hardening buildflags
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} $ENV{CXXFLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} $ENV{LDFLAGS}")

# Creating entries for target: prepair
# ############################

if ( AS_LIBRARY )
  add_library( prepair SHARED TriangleInfo.cpp PolygonRepair.cpp )
else()
  add_executable( prepair TriangleInfo.cpp PolygonRepair.cpp prepair.cpp )
endif()

add_to_cached_list( CGAL_EXECUTABLE_TARGETS prepair )

# Link to CGAL and third-party libraries
target_link_libraries(prepair   ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${GDAL_LIBRARY} ${Boost_LIBRARIES})