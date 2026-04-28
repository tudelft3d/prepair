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

cmake_minimum_required(VERSION 3.5)
project( prepair )

set( CMAKE_BUILD_TYPE "Release")

# CGAL (needs >=5.5 for the Polygon_repair package)
find_package( CGAL 6.0 REQUIRED )

# Boost (for program_options CLI parsing)
find_package( Boost 1.50.0 REQUIRED COMPONENTS program_options filesystem)

# GDAL (for I/O)
find_package( GDAL 3.0 REQUIRED )
include_directories( ${GDAL_INCLUDE_DIR} )

add_executable( prepair prepair.cpp )
set_property(TARGET prepair PROPERTY CXX_STANDARD 17)
target_link_libraries(prepair CGAL::CGAL ${GDAL_LIBRARY} ${Boost_LIBRARIES})
install(TARGETS prepair DESTINATION bin)
