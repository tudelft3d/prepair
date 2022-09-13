# prepair
#
# Copyright © 2009-2022,
# Ken Arroyo Ohori    k.ohori@tudelft.nl
# Hugo Ledoux         h.ledoux@tudelft.nl
# Martijn Meijers     b.m.meijers@tudelft.nl
# All rights reserved.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

project( prepair )

cmake_minimum_required(VERSION 3.1)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

# CGAL
find_package( CGAL REQUIRED )

if ( NOT CGAL_FOUND )
  message(SEND_ERROR "prepair requires the CGAL library")
  return()  
endif()

# include helper file
include( ${CGAL_USE_FILE} )

# Boost
find_package( Boost 1.50.0 REQUIRED COMPONENTS program_options filesystem)

if ( NOT Boost_FOUND )
  message(SEND_ERROR "prepair requires the Boost library version 1.50 or higher")
  return()  
endif()

INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

# GDAL
find_package( GDAL 3.0 )

if ( NOT GDAL_FOUND )
  message(SEND_ERROR "prepair requires the GDAL library")
endif()

include_directories( ${GDAL_INCLUDE_DIR} )

# Creating entries for target: prepair
# ############################

add_executable( prepair prepair.cpp )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS prepair )

# Link to CGAL and third-party libraries
target_link_libraries(prepair ${GDAL_LIBRARY} ${Boost_LIBRARIES})