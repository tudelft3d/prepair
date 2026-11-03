# prepair
#
# Copyright © 2009-2026,
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

cmake_minimum_required(VERSION 3.20)

project(prepair
  VERSION 3.0.1
  DESCRIPTION "Automatic repair of invalid GIS polygons"
  HOMEPAGE_URL "https://github.com/tudelft3d/prepair"
  LANGUAGES CXX)

# Default to Release for single-config generators, but never override a build
# type explicitly chosen by the user (e.g. -DCMAKE_BUILD_TYPE=Debug).
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    Debug Release RelWithDebInfo MinSizeRel)
endif()

# CGAL (needs 6.0 for the Polygon_repair package)
find_package(CGAL 6.0 REQUIRED)

# Boost (for program_options CLI parsing)
# Use CONFIG mode: the FindBoost module is deprecated as of CMake 3.30.
find_package(Boost 1.70 REQUIRED COMPONENTS program_options CONFIG)

# GDAL (for I/O)
find_package(GDAL 3.0 REQUIRED)

add_executable(prepair prepair.cpp)
target_compile_features(prepair PRIVATE cxx_std_17)
set_target_properties(prepair PROPERTIES CXX_EXTENSIONS OFF)

target_link_libraries(prepair PRIVATE
  CGAL::CGAL
  Boost::program_options
  GDAL::GDAL)

if(WIN32)
  target_compile_definitions(prepair PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
endif()
if(MSVC)
  # CGAL templates hit MSVC's default object-file section limit.
  target_compile_options(prepair PRIVATE /bigobj /utf-8)
endif()

include(GNUInstallDirs)
install(TARGETS prepair RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES README.md LICENSE.txt DESTINATION ${CMAKE_INSTALL_DOCDIR})

include(CTest)
if(BUILD_TESTING)
  add_test(NAME repair_bowtie
    COMMAND prepair -w "POLYGON((0 0, 0 10, 10 0, 10 10, 0 0))")
  set_tests_properties(repair_bowtie PROPERTIES
    PASS_REGULAR_EXPRESSION "MULTIPOLYGON"
    FAIL_REGULAR_EXPRESSION "Error")

  add_test(NAME repair_unclosed_ring
    COMMAND prepair -w "POLYGON((0 0, 10 0, 10 10, 0 10))")
  set_tests_properties(repair_unclosed_ring PROPERTIES
    PASS_REGULAR_EXPRESSION "^POLYGON"
    FAIL_REGULAR_EXPRESSION "MULTIPOLYGON|Error")

  add_test(NAME repair_valid_polygon
    COMMAND prepair -w "POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))")
  set_tests_properties(repair_valid_polygon PROPERTIES
    PASS_REGULAR_EXPRESSION "^POLYGON"
    FAIL_REGULAR_EXPRESSION "MULTIPOLYGON|Error")

  add_test(NAME repair_3d
    COMMAND prepair -w "POLYGON((0 0 0, 0 10 5, 10 0 10, 10 10 15, 0 0 0))")
  set_tests_properties(repair_3d PROPERTIES
    PASS_REGULAR_EXPRESSION "MULTIPOLYGON"
    FAIL_REGULAR_EXPRESSION "Error")

  add_test(NAME help
    COMMAND prepair --help)
  set_tests_properties(help PROPERTIES
    PASS_REGULAR_EXPRESSION "usage: prepair")
endif()
