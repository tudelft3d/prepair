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
//#define COORDS_3D

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include "Triangle_info.h"

// STL
#include <fstream>

// OGR
#include <gdal/ogrsf_frmts.h>

// CGAL
#ifdef EXACT_CONSTRUCTIONS
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#else
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#endif
#ifdef COORDS_3D
#include <CGAL/Projection_traits_xy_3.h>
#endif
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

namespace prepair {
  
  // Kernel
#ifdef EXACT_CONSTRUCTIONS
  typedef CGAL::Exact_predicates_exact_constructions_kernel TK;
  typedef CGAL::Exact_intersections_tag IT;
#else
  typedef CGAL::Exact_predicates_inexact_constructions_kernel TK;
  typedef CGAL::Exact_predicates_tag IT;
#endif
  
#ifdef COORDS_3D
  typedef CGAL::Projection_traits_xy_3<TK> K;
#else
  typedef TK K;
#endif
  
  typedef CGAL::Triangulation_vertex_base_2<K> TVB;
  typedef CGAL::Triangulation_hierarchy_vertex_base_2<TVB> VB;
  typedef CGAL::Constrained_triangulation_face_base_2<K> FB;
  typedef CGAL::Triangulation_face_base_with_info_2<Triangle_info, K, FB> FBWI;
  typedef CGAL::Triangulation_data_structure_2<VB, FBWI> TDS;
  typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, IT> CDT;
  typedef CGAL::Triangulation_hierarchy_2<CDT> CDTH;
  typedef CGAL::Constrained_triangulation_plus_2<CDTH> Triangulation;
  
  typedef K::Point_2 Point;
  typedef K::Vector_2 Vector;
}

#endif
