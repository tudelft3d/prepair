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
  static void ogr_to_multi_polygon(OGRGeometry *in_geometry, prepair::Multi_polygon<prepair::Point> &out_geometry);
  static OGRGeometry *multi_polygon_to_ogr(prepair::Multi_polygon<prepair::Point> &in_geometry);
  
  template <class Geometry>
  void repair_odd_even(Geometry &in_geometry, prepair::Multi_polygon<prepair::Point> &out_geometry, bool time_results = false);
  void repair_point_set(prepair::Linear_ring<prepair::Point> &in_geometry, prepair::Multi_polygon<prepair::Point> &out_geometry, bool time_results = false);
  void repair_point_set(prepair::Polygon<prepair::Point> &in_geometry, prepair::Multi_polygon<prepair::Point> &out_geometry, bool time_results = false);
  void repair_point_set(prepair::Multi_polygon<prepair::Point> &in_geometry, prepair::Multi_polygon<prepair::Point> &out_geometry, bool time_results = false);
private:
  prepair::Triangulation triangulation;
  
  void insert_constraints(prepair::Triangulation &triangulation, prepair::Linear_ring<prepair::Point> &in_geometry, bool remove_overlapping_constraints = true);
  void insert_constraints(prepair::Triangulation &triangulation, prepair::Polygon<prepair::Point> &in_geometry, bool remove_overlapping_constraints = true);
  void insert_constraints(prepair::Triangulation &triangulation, prepair::Multi_polygon<prepair::Point> &in_geometry, bool remove_overlapping_constraints = true);
  void tag_odd_even(prepair::Triangulation &triangulation);
  void tag_point_set(prepair::Triangulation &triangulation, std::list<prepair::Multi_polygon<prepair::Point> > &geometries, std::list<bool> &geometries_flipped);
  void reconstruct(prepair::Triangulation &triangulation, prepair::Multi_polygon<prepair::Point> &out_geometry);
  void get_boundary(prepair::Triangulation::Face_handle face, int edge, prepair::Linear_ring<prepair::Triangulation::Vertex_handle> &out_vertices);
  
  // Debug functions
  void print_triangle(prepair::Triangulation::Face_handle triangle);
};

// NOTE: Has to be here because C++'s compilation model doesn't allow templated classes in implementation files
template <class Geometry>
void Polygon_repair::repair_odd_even(Geometry &in_geometry, prepair::Multi_polygon<prepair::Point> &out_geometry, bool time_results) {
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