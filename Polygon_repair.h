/*
 Copyright (c) 2009-2014,
 Ken Arroyo Ohori    g.a.k.arroyoohori@tudelft.nl
 Hugo Ledoux         h.ledoux@tudelft.nl
 Martijn Meijers     b.m.meijers@tudelft.nl
 All rights reserved.
 
 This file is part of prepair: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Licensees holding a valid commercial license may use this file in
 accordance with the commercial license agreement provided with
 the software.
 
 This file is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef POLYGONREPAIR_H
#define POLYGONREPAIR_H

#include "definitions.h"

class Polygon_repair {
public:
  static void ogr_to_multi_polygon(OGRGeometry *in_geometry, Multi_polygon<Point_2> &out_geometry);
  static OGRGeometry *multi_polygon_to_ogr(Multi_polygon<Point_2> &in_geometry);
  
  template <class Geometry>
  void repair_odd_even(Geometry &in_geometry, Multi_polygon<Point_2> &out_geometry, bool time_results = false);
  void repair_point_set(Linear_ring<Point_2> &in_geometry, Multi_polygon<Point_2> &out_geometry, bool time_results = false);
  void repair_point_set(Polygon<Point_2> &in_geometry, Multi_polygon<Point_2> &out_geometry, bool time_results = false);
  void repair_point_set(Multi_polygon<Point_2> &in_geometry, Multi_polygon<Point_2> &out_geometry, bool time_results = false);
private:
  Triangulation triangulation;
  
  void insert_constraints(Triangulation &triangulation, Linear_ring<Point_2> &in_geometry, bool remove_overlapping_constraints = true);
  void insert_constraints(Triangulation &triangulation, Polygon<Point_2> &in_geometry, bool remove_overlapping_constraints = true);
  void insert_constraints(Triangulation &triangulation, Multi_polygon<Point_2> &in_geometry, bool remove_overlapping_constraints = true);
  void tag_odd_even(Triangulation &triangulation);
  void tag_point_set(Triangulation &triangulation, std::list<Multi_polygon<Point_2> > &geometries, std::list<bool> &geometries_flipped);
  void reconstruct(Triangulation &triangulation, Multi_polygon<Point_2> &out_geometry);
  void get_boundary(Triangulation::Face_handle face, int edge, Linear_ring<Triangulation::Vertex_handle> &out_vertices);
  
  // Debug functions
  void print_triangle(Triangulation::Face_handle triangle);
};

// NOTE: Has to be here because C++'s compilation model doesn't allow templated classes in implementation files
template <class Geometry>
void Polygon_repair::repair_odd_even(Geometry &in_geometry, Multi_polygon<Point_2> &out_geometry, bool time_results) {
  triangulation.clear();
  time_t this_time, total_time;
  this_time = time(NULL);
  insert_constraints(triangulation, in_geometry);
  total_time = time(NULL)-this_time;
  if (time_results) std::cout << "Triangulation: " << total_time/60 << " minutes " << total_time%60 << " seconds." << std::endl;
  this_time = time(NULL);
  tag_odd_even(triangulation);
  total_time = time(NULL)-this_time;
  if (time_results) std::cout << "Tagging: " << total_time/60 << " minutes " << total_time%60 << " seconds." << std::endl;
  this_time = time(NULL);
  reconstruct(triangulation, out_geometry);
  total_time = time(NULL)-this_time;
  if (time_results) std::cout << "Reconstruction: " << total_time/60 << " minutes " << total_time%60 << " seconds." << std::endl;
}

#endif