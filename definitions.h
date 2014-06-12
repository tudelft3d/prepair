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

// Compile-time options
// * if the code crashes, try to compile with EXACT_CONSTRUCTIONS so that
//   robust arithmetic is used

#define EXACT_CONSTRUCTIONS
#define COORDS_3D

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include "Triangle_info.h"

// STL
#include <iostream>
#include <stack>
#include <set>
#include <list>
#include <fstream>
#include <math.h>
#include <time.h>

// OGR
#include <gdal/ogrsf_frmts.h>

// CGAL
#ifdef EXACT_CONSTRUCTIONS
  #include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#else
  #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#endif
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Snap_rounding_traits_2.h>
#include <CGAL/Snap_rounding_2.h>

// Kernel
#ifdef EXACT_CONSTRUCTIONS
  typedef CGAL::Exact_predicates_exact_constructions_kernel K;
#else
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
#endif

typedef CGAL::Triangulation_vertex_base_2<K> TVB;
typedef CGAL::Triangulation_hierarchy_vertex_base_2<TVB> VB;
typedef CGAL::Constrained_triangulation_face_base_2<K> FB;
typedef CGAL::Triangulation_face_base_with_info_2<Triangle_info, K, FB> FBWI;
typedef CGAL::Triangulation_data_structure_2<VB, FBWI> TDS;
typedef CGAL::Exact_predicates_tag PT;
typedef CGAL::Exact_intersections_tag IT;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, PT> CDT;
typedef CGAL::Triangulation_hierarchy_2<CDT> CDTH;
typedef CGAL::Constrained_triangulation_plus_2<CDTH> Triangulation;
typedef CGAL::Snap_rounding_traits_2<K> SRT;

typedef K::Point_2 Point;
typedef K::Segment_2 Segment;
typedef K::Vector_2 Vector;

// Custom containers

template <class Vertex>
class Vertex_converter;

template <class Vertex_>
class Linear_ring {
private:
  std::list<Vertex_> v;
public:
  typedef Vertex_ Vertex;
  typedef typename std::list<Vertex>::iterator Vertex_iterator;
  typedef typename std::list<Vertex>::reverse_iterator Vertex_reverse_iterator;
  
  void clear() {
    v.clear();
  }
  
  std::size_t number_of_vertices() {
    return v.size();
  }
  
  Vertex &first_vertex() {
    return v.front();
  }
  
  Vertex &last_vertex() {
    return v.back();
  }
  
  void add_vertex(Vertex vertex) {
    v.push_back(vertex);
  }
  
  Vertex_iterator vertices_begin() {
    return v.begin();
  }
  
  Vertex_iterator vertices_end() {
    return v.end();
  }
  
  Vertex_reverse_iterator vertices_rbegin() {
    return v.rbegin();
  }
  
  Vertex_reverse_iterator vertices_rend() {
    return v.rend();
  }
  
  void splice_vertices(Vertex_iterator position, Linear_ring &other) {
    v.splice(position, other.v);
  }
  
  bool is_clockwise() {
    Vertex_iterator current = vertices_begin();
    Vertex_iterator rightmost = current;
    while (current != vertices_end()) {
      if (Vertex_converter<Vertex>::to_point(*current).x() > Vertex_converter<Vertex>::to_point(*rightmost).x()) {
        rightmost = current;
      } ++current;
    } Vertex_iterator next = rightmost;
    ++next;
    if (next == vertices_end()) next = vertices_begin();
    Vertex_iterator previous = rightmost;
    if (previous == vertices_begin()) previous = vertices_end();
    --previous;
    
    if (Vertex_converter<Vertex>::to_point(*previous).y() < Vertex_converter<Vertex>::to_point(*rightmost).y() &&
        Vertex_converter<Vertex>::to_point(*rightmost).y() < Vertex_converter<Vertex>::to_point(*next).y()) return true;
    
    return false;
  }
};

template <class Vertex_>
class Polygon {
private:
  Linear_ring<Vertex_> oring;
  std::list<Linear_ring<Vertex_> > irings;
public:
  typedef Vertex_ Vertex;
  typedef typename std::list<Linear_ring<Vertex> >::iterator Inner_ring_iterator;
  typedef typename std::list<Linear_ring<Vertex> >::reverse_iterator Inner_ring_reverse_iterator;
  
  void clear() {
    oring.clear();
    irings.clear();
  }
  
  std::size_t number_of_inner_rings() {
    return irings.size();
  }
  
  Linear_ring<Vertex> &outer_ring() {
    return oring;
  }
  
  Linear_ring<Vertex> &first_inner_ring() {
    return irings.front();
  }
  
  Linear_ring<Vertex> &last_inner_ring() {
    return irings.back();
  }
  
  Linear_ring<Vertex> &add_empty_inner_ring() {
    irings.push_back(Linear_ring<Vertex>());
    return irings.back();
  }
  
  Inner_ring_iterator inner_rings_begin() {
    return irings.begin();
  }
  
  Inner_ring_iterator inner_rings_end() {
    return irings.end();
  }
  
  Inner_ring_reverse_iterator inner_rings_rbegin() {
    return irings.rbegin();
  }
  
  Inner_ring_reverse_iterator inner_rings_rend() {
    return irings.rend();
  }
};

template <class Vertex_>
class Multi_polygon {
private:
  std::list<Polygon<Vertex_> > p;
public:
  typedef Vertex_ Vertex;
  typedef typename std::list<Polygon<Vertex> >::iterator Polygon_iterator;
  typedef typename std::list<Polygon<Vertex> >::reverse_iterator Polygon_reverse_iterator;
  
  void clear() {
    p.clear();
  }
  
  std::size_t number_of_polygons() {
    return p.size();
  }
  
  Polygon<Vertex> &first_polygon() {
    return p.front();
  }
  
  Polygon<Vertex> &last_polygon() {
    return p.back();
  }
  
  Polygon<Vertex> &add_empty_polygon() {
    p.push_back(Polygon<Vertex>());
    return p.back();
  }
  
  Polygon_iterator polygons_begin() {
    return p.begin();
  }
  
  Polygon_iterator polygons_end() {
    return p.end();
  }
  
  Polygon_reverse_iterator polygons_rbegin() {
    return p.rbegin();
  }
  
  Polygon_reverse_iterator polygons_rend() {
    return p.rend();
  }
  
  std::string as_wkt() {
    std::string wkt = "MULTIPOLYGON(";
    for (Polygon_iterator current_polygon = polygons_begin(); current_polygon != polygons_end(); ++current_polygon) {
      wkt.append("((");
      for (typename Linear_ring<Vertex>::Vertex_iterator current_vertex = current_polygon->outer_ring().vertices_begin(); current_vertex != current_polygon->outer_ring().vertices_end(); ++current_vertex) {
        if (current_vertex != current_polygon->outer_ring().vertices_begin()) wkt.append(", ");
        wkt.append(std::to_string(CGAL::to_double(Vertex_converter<Vertex>::to_point(*current_vertex).x())));
        wkt.append(" ");
        wkt.append(std::to_string(CGAL::to_double(Vertex_converter<Vertex>::to_point(*current_vertex).y())));
      } for (typename Polygon<Vertex>::Inner_ring_iterator current_ring = current_polygon->inner_rings_begin(); current_ring != current_polygon->inner_rings_end(); ++current_ring) {
        if (current_ring != current_polygon->inner_rings_begin()) wkt.append(", ");
        for (typename Linear_ring<Vertex>::Vertex_iterator current_vertex = current_ring->vertices_begin(); current_vertex != current_ring->vertices_end(); ++current_vertex) {
          if (current_vertex != current_ring->vertices_begin()) wkt.append(", ");
          wkt.append(std::to_string(CGAL::to_double(Vertex_converter<Vertex>::to_point(*current_vertex).x())));
          wkt.append(" ");
          wkt.append(std::to_string(CGAL::to_double(Vertex_converter<Vertex>::to_point(*current_vertex).y())));
        }
      } wkt.append(")");
    } wkt.append(")");
    return wkt;
  }
};

template <>
class Vertex_converter<Triangulation::Vertex_handle> {
public:
  static Point to_point(Triangulation::Vertex_handle vertex) {
    return vertex->point();
  }
};

template <>
class Vertex_converter<Point> {
public:
  static Point to_point(Point vertex) {
    return vertex;
  }
};

#endif
