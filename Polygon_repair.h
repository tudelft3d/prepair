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
  typedef prepair::Triangulation Triangulation;
  typedef prepair::Point Point;
  typedef prepair::Vector Vector;
  
  bool is_iso_and_ogc_valid(OGRGeometry *in_geometry, const std::string &pre_text = std::string("\t"), bool time_results = false);
  OGRGeometry *repair_odd_even(OGRGeometry *in_geometry, bool time_results = false);
  OGRGeometry *repair_point_set(OGRGeometry *in_geometry, bool time_results = false);
  void remove_small_parts(OGRGeometry *geometry, double min_area);

private:
  Triangulation triangulation;
  Triangulation::Face_handle walk_start_location;
  
  void insert_constraints(OGRGeometry *in_geometry);
  void tag_odd_even();
  void tag_point_set_difference(std::list<std::pair<bool, OGRGeometry *> > &geometries);
  void tag_point_set_union(std::list<std::pair<bool, OGRGeometry *> > &geometries);
  OGRGeometry *reconstruct();
  void get_boundary(Triangulation::Face_handle face, int edge, std::list<Triangulation::Vertex_handle> &out_vertices);
};

#endif