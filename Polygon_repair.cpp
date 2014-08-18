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

#include "Polygon_repair.h"

bool Polygon_repair::is_iso_and_ogc_valid(OGRGeometry *in_geometry, const std::string &pre_text, bool time_results) {
  // TODO: Implement
  triangulation.clear();
  std::time_t this_time, total_time;
  bool is_valid = true;
  
  this_time = time(NULL);
  switch (in_geometry->getGeometryType()) {
    case wkbLineString: {
      OGRLinearRing *ring = static_cast<OGRLinearRing *>(in_geometry);
      
      if (ring->IsEmpty()) {
        std::cout << pre_text << "ring: empty" << std::endl;
        is_valid = false;
      }
      
      if (is_valid) {
#ifdef COORDS_3D
        Point first(ring->getX(0), ring->getY(0), ring->getZ(0));
        Point last(ring->getX(ring->getNumPoints()-1), ring->getY(ring->getNumPoints()-1), ring->getZ(ring->getNumPoints()-1));
#else
        Point first(ring->getX(0), ring->getY(0));
        Point last(ring->getX(ring->getNumPoints()-1), ring->getY(ring->getNumPoints()-1));
#endif
        if (first != last) {
          std::cout << pre_text << "ring: not closing edge between (" << last << " and " << first << ")" << std::endl;
          is_valid = false;
          ring->closeRings();
        }
      }
      
      if (is_valid) {
        if (ring->getNumPoints() < 4) {
          std::cout << pre_text << "ring: less than 4 vertices" << std::endl;
          is_valid = false;
        }
      }
      
      if (is_valid) {
#ifdef COORDS_3D
        Point previous(ring->getX(0), ring->getY(0), ring->getZ(0));
#else
        Point previous(ring->getX(0), ring->getY(0));
#endif
        for (int currentVertex = 1; currentVertex < ring->getNumPoints(); ++currentVertex) {
#ifdef COORDS_3D
          Point current(ring->getX(currentVertex), ring->getY(currentVertex), ring->getZ(currentVertex));
#else
          Point current(ring->getX(currentVertex), ring->getY(currentVertex));
#endif
          if (previous == current) {
            std::cout << pre_text << "ring vertex " << currentVertex << ": duplicate point" << std::endl;
          }
        }
      }
      
      total_time = time(NULL)-this_time;
      if (time_results) std::cout << "Simple checks: " << total_time/60 << " minutes " << total_time%60 << " seconds." << std::endl;
      if (!is_valid) return false;
      
      this_time = time(NULL);
      
      total_time = time(NULL)-this_time;
      if (time_results) std::cout << "Triangulation checks: " << total_time/60 << " minutes " << total_time%60 << " seconds." << std::endl;
      
      
    }
      
      
      
      
      
      
    
      
    default:
      std::cerr << "Error: Input type not supported" << std::endl;
      return false;
      break;
  }
  
  
  
}

OGRGeometry *Polygon_repair::repair_odd_even(OGRGeometry *in_geometry, bool time_results) {
  triangulation.clear();
  std::time_t this_time, total_time;
  
  this_time = time(NULL);
  insert_odd_even_constraints(in_geometry);
  total_time = time(NULL)-this_time;
  if (time_results) std::cout << "Triangulation: " << total_time/60 << " minutes " << total_time%60 << " seconds." << std::endl;
  
  this_time = time(NULL);
  tag_odd_even();
  total_time = time(NULL)-this_time;
  if (time_results) std::cout << "Tagging: " << total_time/60 << " minutes " << total_time%60 << " seconds." << std::endl;
  
  this_time = time(NULL);
  OGRGeometry *out_geometry = reconstruct();
  total_time = time(NULL)-this_time;
  if (time_results) std::cout << "Reconstruction: " << total_time/60 << " minutes " << total_time%60 << " seconds." << std::endl;
  
  return out_geometry;
}

OGRGeometry *Polygon_repair::repair_point_set(OGRGeometry *in_geometry, bool time_results) {
  std::time_t this_time, total_time;
  std::list<OGRGeometry *> repaired_parts;   // bool indicates if outer/inner are flipped
  
  this_time = time(NULL);
  switch (in_geometry->getGeometryType()) {
    case wkbLineString: {
      return repair_odd_even(in_geometry, time_results);
      break;
    }
      
    case wkbPolygon: {
      OGRPolygon *polygon = static_cast<OGRPolygon *>(in_geometry);
      repaired_parts.push_back(repair_point_set(polygon->getExteriorRing()));
      for (int current_ring = 0; current_ring < polygon->getNumInteriorRings(); ++current_ring) {
        repaired_parts.push_back(repair_point_set(polygon->getInteriorRing(current_ring)));
      } total_time = time(NULL)-this_time;
      if (time_results) std::cout << "Repairing individual rings: " << total_time/60 << " minutes " << total_time%60 << " seconds." << std::endl;
      
      this_time = time(NULL);
      triangulation.clear();
      for (std::list<OGRGeometry *>::iterator current_part = repaired_parts.begin(); current_part != repaired_parts.end(); ++current_part) {
        insert_all_constraints(*current_part);
      } total_time = time(NULL)-this_time;
      if (time_results) std::cout << "Triangulation: " << total_time/60 << " minutes " << total_time%60 << " seconds." << std::endl;
      
      this_time = time(NULL);
      tag_point_set_difference(repaired_parts);
      total_time = time(NULL)-this_time;
      if (time_results) std::cout << "Tagging: " << total_time/60 << " minutes " << total_time%60 << " seconds." << std::endl;
      
      break;
    }
      
    case wkbMultiPolygon: {
      OGRMultiPolygon *multipolygon = static_cast<OGRMultiPolygon *>(in_geometry);
      for (int current_polygon = 0; current_polygon < multipolygon->getNumGeometries(); ++current_polygon) {
        OGRPolygon *polygon = static_cast<OGRPolygon *>(multipolygon->getGeometryRef(current_polygon));
        repaired_parts.push_back(repair_point_set(polygon));
      } total_time = time(NULL)-this_time;
      if (time_results) std::cout << "Repairing individual polygons: " << total_time/60 << " minutes " << total_time%60 << " seconds." << std::endl;
      
      this_time = time(NULL);
      triangulation.clear();
      for (std::list<OGRGeometry *>::iterator current_part = repaired_parts.begin(); current_part != repaired_parts.end(); ++current_part) {
        insert_all_constraints(*current_part);
      } total_time = time(NULL)-this_time;
      if (time_results) std::cout << "Triangulation: " << total_time/60 << " minutes " << total_time%60 << " seconds." << std::endl;
      
      this_time = time(NULL);
      tag_point_set_union(repaired_parts);
      total_time = time(NULL)-this_time;
      if (time_results) std::cout << "Tagging: " << total_time/60 << " minutes " << total_time%60 << " seconds." << std::endl;
      
      break;
    }
      
    default:
      std::cerr << "Error: Input type not supported" << std::endl;
      return new OGRPolygon();
      break;
  }
  
  this_time = time(NULL);
  OGRGeometry *out_geometry = reconstruct();
  total_time = time(NULL)-this_time;
  if (time_results) std::cout << "Reconstruction: " << total_time/60 << " minutes " << total_time%60 << " seconds." << std::endl;
  
  return out_geometry;
}

void Polygon_repair::insert_all_constraints(OGRGeometry *in_geometry) {
  Triangulation::Vertex_handle va, vb;
  
  switch (in_geometry->getGeometryType()) {
    case wkbLineString: {
      OGRLinearRing *ring = static_cast<OGRLinearRing *>(in_geometry);
      ring->closeRings();
      
#ifdef COORDS_3D
      vb = triangulation.insert(Point(ring->getX(0), ring->getY(0), ring->getZ(0)), walk_start_location);
#else
      vb = triangulation.insert(Point(ring->getX(0), ring->getY(0)), walk_start_location);
#endif
      walk_start_location = triangulation.incident_faces(vb);
      for (int current_point = 1; current_point < ring->getNumPoints(); ++current_point) {
        va = vb;
#ifdef COORDS_3D
        vb = triangulation.insert(Point(ring->getX(current_point),
                                        ring->getY(current_point),
                                        ring->getZ(current_point)),
                                  walk_start_location);
#else
        vb = triangulation.insert(Point(ring->getX(current_point),
                                        ring->getY(current_point)),
                                  walk_start_location);
#endif
        if (va == vb) continue;
        triangulation.insert_constraint(va, vb);
        walk_start_location = triangulation.incident_faces(vb);
      } break;
    }
      
    case wkbPolygon: {
      OGRPolygon *polygon = static_cast<OGRPolygon *>(in_geometry);
      insert_odd_even_constraints(polygon->getExteriorRing());
      for (int current_ring = 0; current_ring < polygon->getNumInteriorRings(); ++current_ring) {
        insert_all_constraints(polygon->getInteriorRing(current_ring));
      } break;
    }
      
    case wkbMultiPolygon: {
      OGRMultiPolygon *multipolygon = static_cast<OGRMultiPolygon *>(in_geometry);
      for (int current_polygon = 0; current_polygon < multipolygon->getNumGeometries(); ++current_polygon) {
        insert_all_constraints(multipolygon->getGeometryRef(current_polygon));
      } break;
    }
      
    default:
      std::cerr << "Error: Input type not supported" << std::endl;
      return;
      break;
  }
}

void Polygon_repair::insert_odd_even_constraints(OGRGeometry *in_geometry) {
  Triangulation::Vertex_handle va, vb;
  
  switch (in_geometry->getGeometryType()) {
    case wkbLineString: {
      OGRLinearRing *ring = static_cast<OGRLinearRing *>(in_geometry);
      ring->closeRings();
#ifdef COORDS_3D
      vb = triangulation.insert(Point(ring->getX(0), ring->getY(0), ring->getZ(0)), walk_start_location);
#else
      vb = triangulation.insert(Point(ring->getX(0), ring->getY(0)), walk_start_location);
#endif
      walk_start_location = triangulation.incident_faces(vb);
      for (int current_point = 1; current_point < ring->getNumPoints(); ++current_point) {
        va = vb;
#ifdef COORDS_3D
        vb = triangulation.insert(Point(ring->getX(current_point),
                                        ring->getY(current_point),
                                        ring->getZ(current_point)),
                                  walk_start_location);
#else
        vb = triangulation.insert(Point(ring->getX(current_point),
                                        ring->getY(current_point)),
                                  walk_start_location);
#endif
        if (va == vb) continue;
        triangulation.odd_even_insert_constraint(va, vb);
        walk_start_location = triangulation.incident_faces(vb);
      } break;
    }
      
    case wkbPolygon: {
      OGRPolygon *polygon = static_cast<OGRPolygon *>(in_geometry);
      insert_odd_even_constraints(polygon->getExteriorRing());
      for (int current_ring = 0; current_ring < polygon->getNumInteriorRings(); ++current_ring) {
        insert_odd_even_constraints(polygon->getInteriorRing(current_ring));
      } break;
    }
      
    case wkbMultiPolygon: {
      OGRMultiPolygon *multipolygon = static_cast<OGRMultiPolygon *>(in_geometry);
      for (int current_polygon = 0; current_polygon < multipolygon->getNumGeometries(); ++current_polygon) {
        insert_odd_even_constraints(multipolygon->getGeometryRef(current_polygon));
      } break;
    }
      
    default:
      std::cerr << "Error: Input type not supported" << std::endl;
      return;
      break;
  }
}

void Polygon_repair::tag_odd_even() {
  // Clean tags
  for (prepair::Triangulation::Face_handle current_face = triangulation.all_faces_begin(); current_face != triangulation.all_faces_end(); ++current_face)
    current_face->info().clear();
  
  // Initialise tagging
  std::stack<prepair::Triangulation::Face_handle> interior_stack, exterior_stack;
  exterior_stack.push(triangulation.infinite_face());
  std::stack<prepair::Triangulation::Face_handle> *current_stack = &exterior_stack;
  std::stack<prepair::Triangulation::Face_handle> *dual_stack = &interior_stack;
  bool labelling_interior = false;
  
  
  // Until we finish
  while (!interior_stack.empty() || !exterior_stack.empty()) {
    
    // Give preference to whatever we're already doing
    while (!current_stack->empty()) {
      prepair::Triangulation::Face_handle current_face = current_stack->top();
			current_stack->pop();
      if (current_face->info().been_tagged()) continue;
			current_face->info().is_in_interior(labelling_interior);
      for (int current_edge = 0; current_edge < 3; ++current_edge) {
        if (!current_face->neighbor(current_edge)->info().been_tagged()) {
          if (current_face->is_constrained(current_edge))
            dual_stack->push(current_face->neighbor(current_edge));
          else
            current_stack->push(current_face->neighbor(current_edge));
        }
      }
    }
    
    // Flip
    if (!labelling_interior) {
      current_stack = &interior_stack;
      dual_stack = &exterior_stack;
    } else {
      current_stack = &exterior_stack;
      dual_stack = &interior_stack;
    } labelling_interior = !labelling_interior;
	}
}

void Polygon_repair::tag_as_to_fill_in(OGRGeometry *geometry) {
  Triangulation::Vertex_handle va, vb;
  Triangulation::Face_handle face;
  int index_of_opposite_vertex;
  
  switch (geometry->getGeometryType()) {
      
    case wkbLineString: {
      OGRLinearRing *ring = static_cast<OGRLinearRing *>(geometry);
#ifdef COORDS_3D
      vb = triangulation.insert(Point(ring->getX(0), ring->getY(0), ring->getZ(0)), walk_start_location);
#else
      vb = triangulation.insert(Point(ring->getX(0), ring->getY(0)), walk_start_location);
#endif
      walk_start_location = triangulation.incident_faces(vb);
      for (int current_point = 1; current_point < ring->getNumPoints(); ++current_point) {
        va = vb;
#ifdef COORDS_3D
        vb = triangulation.insert(Point(ring->getX(current_point),
                                        ring->getY(current_point),
                                        ring->getZ(current_point)),
                                  walk_start_location);
#else
        vb = triangulation.insert(Point(ring->getX(current_point),
                                        ring->getY(current_point)),
                                  walk_start_location);
#endif
        if (triangulation.is_edge(va, vb, face, index_of_opposite_vertex)) {
          
        }
        walk_start_location = triangulation.incident_faces(vb);
      } break;
    }
      
    default:
      std::cerr << "Error: Input type not supported" << std::endl;
      return;
      break;
  }
}

void Polygon_repair::tag_as_to_carve_out(OGRGeometry *geometry) {
  
}

void Polygon_repair::tag_point_set_difference(std::list<OGRGeometry *> &geometries) {
  // TODO: Implement
  std::list<OGRGeometry *>::iterator current_geometry = geometries.begin();
  ++current_geometry;
  tag_as_to_fill_in(*current_geometry);
  while (current_geometry != geometries.end()) {
    tag_as_to_carve_out(*current_geometry);
    ++current_geometry;
  }
}

void Polygon_repair::tag_point_set_union(std::list<OGRGeometry *> &geometries) {
  // TODO: Implement
}

OGRGeometry *Polygon_repair::reconstruct() {
  // std::cout << "Triangulation: " << triangulation.number_of_faces() << " faces, " << triangulation.number_of_vertices() << " vertices." << std::endl;
  if (triangulation.number_of_faces() < 1) {
    return new OGRPolygon();
  }
  
  // Reconstruct
  OGRMultiPolygon *out_geometries = new OGRMultiPolygon();
  for (prepair::Triangulation::Finite_faces_iterator seeding_face = triangulation.finite_faces_begin(); seeding_face != triangulation.finite_faces_end(); ++seeding_face) {
    
    if (!seeding_face->info().is_in_interior() || seeding_face->info().been_reconstructed()) continue;
    seeding_face->info().been_reconstructed(true);
    if (!seeding_face->info().been_reconstructed()) {
      std::cout << "Error: Face should be marked as reconstructed!" << std::endl;
    }
    
    // Get boundary
    std::list<Triangulation::Vertex_handle> vertices = std::list<Triangulation::Vertex_handle>();
    if (seeding_face->neighbor(2)->info().is_in_interior() && !seeding_face->neighbor(2)->info().been_reconstructed()) {
      seeding_face->neighbor(2)->info().been_reconstructed(true);
      std::list<Triangulation::Vertex_handle> l2;
      get_boundary(seeding_face->neighbor(2), seeding_face->neighbor(2)->index(seeding_face), l2);
      vertices.splice(vertices.end(), l2);
    } vertices.push_back(seeding_face->vertex(0));
    if (seeding_face->neighbor(1)->info().is_in_interior() && !seeding_face->neighbor(1)->info().been_reconstructed()) {
      seeding_face->neighbor(1)->info().been_reconstructed(true);
      std::list<Triangulation::Vertex_handle> l1;
      get_boundary(seeding_face->neighbor(1), seeding_face->neighbor(1)->index(seeding_face), l1);
      vertices.splice(vertices.end(), l1);
    } vertices.push_back(seeding_face->vertex(2));
    if (seeding_face->neighbor(0)->info().is_in_interior() && !seeding_face->neighbor(0)->info().been_reconstructed()) {
      seeding_face->neighbor(0)->info().been_reconstructed(true);
      std::list<Triangulation::Vertex_handle> l0;
      get_boundary(seeding_face->neighbor(0), seeding_face->neighbor(0)->index(seeding_face), l0);
      vertices.splice(vertices.end(), l0);
    } vertices.push_back(seeding_face->vertex(1));
    
    // Find cutting vertices
    std::set<Triangulation::Vertex_handle> visited_vertices;
    std::set<Triangulation::Vertex_handle> repeated_vertices;
    for (std::list<Triangulation::Vertex_handle>::iterator current_vertex = vertices.begin(); current_vertex != vertices.end(); ++current_vertex) {
      if (!visited_vertices.insert(*current_vertex).second) repeated_vertices.insert(*current_vertex);
    } visited_vertices.clear();
    
    // Cut and join rings in the correct order
    std::list<std::list<Triangulation::Vertex_handle> > rings;
    std::stack<std::list<Triangulation::Vertex_handle> > chains_stack;
    std::set<Triangulation::Vertex_handle> vertices_where_chains_begin;
    rings.push_back(std::list<Triangulation::Vertex_handle>());
    for (std::list<Triangulation::Vertex_handle>::iterator current_vertex = vertices.begin(); current_vertex != vertices.end(); ++current_vertex) {
      
      // New chain
      if (repeated_vertices.count(*current_vertex) > 0) {
        // Closed by itself
        if (rings.back().front() == *current_vertex) {
          // Degenerate (insufficient vertices to be valid)
          if (rings.back().size() < 3) {
            rings.back().clear();
          } else {
            std::list<Triangulation::Vertex_handle>::iterator second_element = rings.back().begin();
            ++second_element;
            // Degenerate (zero area)
            if (rings.back().back() == *second_element) {
              rings.back().clear();
            }
            // Valid
            else {
              rings.push_back(std::list<Triangulation::Vertex_handle>());
            }
          }
        }
        // Open by itself
        else {
          // Closed with others in stack
          if (vertices_where_chains_begin.count(*current_vertex)) {
            
            while (rings.back().front() != *current_vertex) {
              rings.back().splice(rings.back().begin(), chains_stack.top());
              chains_stack.pop();
            } vertices_where_chains_begin.erase(*current_vertex);
            // Degenerate (insufficient vertices to be valid)
            if (rings.back().size() < 3) {
              rings.back().clear();
            } else {
              std::list<Triangulation::Vertex_handle>::iterator second_element = rings.back().begin();
              ++second_element;
              // Degenerate (zero area)
              if (rings.back().back() == *second_element) {
                rings.back().clear();
              }
              // Valid
              else {
                rings.push_back(std::list<Triangulation::Vertex_handle>());
              }
            }
          }
          // Open
          else {
            // Not first chain
            if (repeated_vertices.count(rings.back().front()) > 0) {
              vertices_where_chains_begin.insert(rings.back().front());
            }
            chains_stack.push(std::list<Triangulation::Vertex_handle>());
            chains_stack.top().splice(chains_stack.top().begin(), rings.back());
          }
        }
      } rings.back().push_back(*current_vertex);
    }
    // Final ring
    while (chains_stack.size() > 0) {
      rings.back().splice(rings.back().begin(), chains_stack.top());
      chains_stack.pop();
    }
    // Degenerate (insufficient vertices to be valid)
    if (rings.back().size() < 3) {
      rings.back().clear();
    } else {
      std::list<Triangulation::Vertex_handle>::iterator second_element = rings.back().begin();
      ++second_element;
      // Degenerate (zero area)
      if (rings.back().back() == *second_element) {
        rings.back().clear();
      }
    }
    
    // Remove last ring if too small (or empty)
    if (rings.back().size() < 3) {
      rings.pop_back();
    }
    
    // Start rings at the lexicographically smallest vertex
    for (std::list<std::list<Triangulation::Vertex_handle> >::iterator current_ring = rings.begin(); current_ring != rings.end(); ++current_ring) {
      std::list<Triangulation::Vertex_handle>::iterator smallest_vertex = current_ring->begin();
      for (std::list<Triangulation::Vertex_handle>::iterator current_vertex = current_ring->begin(); current_vertex != current_ring->end(); ++current_vertex) {
        if ((*current_vertex)->point() < (*smallest_vertex)->point()) smallest_vertex = current_vertex;
      } if (current_ring->back() != *smallest_vertex) {
        ++smallest_vertex;
        current_ring->splice(current_ring->begin(), *current_ring, smallest_vertex, current_ring->end());
      }
    }
    
    // Make rings
    if (rings.size() == 0) continue;
    std::list<OGRLinearRing *> rings_for_polygon;
    for (std::list<std::list<Triangulation::Vertex_handle> >::iterator current_ring = rings.begin(); current_ring != rings.end(); ++current_ring) {
      OGRLinearRing *new_ring = new OGRLinearRing();
      for (std::list<Triangulation::Vertex_handle>::reverse_iterator current_vertex = current_ring->rbegin(); current_vertex != current_ring->rend(); ++current_vertex) {
        new_ring->addPoint(CGAL::to_double((*current_vertex)->point().x()), CGAL::to_double((*current_vertex)->point().y()));
      } new_ring->addPoint(CGAL::to_double(current_ring->back()->point().x()), CGAL::to_double(current_ring->back()->point().y()));
      rings_for_polygon.push_back(new_ring);
    } OGRPolygon *new_polygon = new OGRPolygon();
    for (std::list<OGRLinearRing *>::iterator current_ring = rings_for_polygon.begin(); current_ring != rings_for_polygon.end(); ++current_ring) {
      if (!(*current_ring)->isClockwise()) {
        new_polygon->addRingDirectly(*current_ring);
        break;
      }
    } for (std::list<OGRLinearRing *>::iterator current_ring = rings_for_polygon.begin(); current_ring != rings_for_polygon.end(); ++current_ring)
      if ((*current_ring)->isClockwise()) new_polygon->addRingDirectly(*current_ring);
    out_geometries->addGeometryDirectly(new_polygon);
  }
  
  if (out_geometries->getNumGeometries() == 0) {
    delete out_geometries;
    return new OGRPolygon();
  }
  
  if (out_geometries->getNumGeometries() == 1) {
    OGRPolygon *new_polygon = static_cast<OGRPolygon *>(out_geometries->getGeometryRef(0)->clone());
    delete out_geometries;
    return new_polygon;
  }
  
  return out_geometries;
}

void remove_small_parts(OGRGeometry *geometry, double min_area) {
  switch (geometry->getGeometryType()) {
    case wkbPolygon: {
      OGRPolygon *polygon = static_cast<OGRPolygon *>(geometry);
      if (polygon->getExteriorRing()->get_Area() < min_area) {
        polygon->empty();
        return;
      } OGRPolygon *new_polygon = new OGRPolygon();
      new_polygon->addRing(polygon->getExteriorRing());
      for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing) {
        if (polygon->getInteriorRing(currentRing)->get_Area() >= min_area) {
          new_polygon->addRing(polygon->getInteriorRing(currentRing));
        }
      } delete geometry;
      geometry = new_polygon;
      break;
    }
      
    case wkbMultiPolygon: {
      OGRMultiPolygon *multipolygon = static_cast<OGRMultiPolygon *>(geometry);
      for (int current_polygon = multipolygon->getNumGeometries()-1; current_polygon > 0; --current_polygon) {
        remove_small_parts(multipolygon->getGeometryRef(current_polygon), current_polygon);
        if (static_cast<OGRPolygon *>(multipolygon->getGeometryRef(current_polygon))->get_Area() < min_area) {
          multipolygon->removeGeometry(current_polygon);
        }
      } break;
    }
      
    default:
      std::cerr << "Error: Input type not supported" << std::endl;
      return;
      break;
  }
}

void Polygon_repair::get_boundary(Triangulation::Face_handle face, int edge, std::list<Triangulation::Vertex_handle> &out_vertices) {
  // Check clockwise edge
  if (face->neighbor(face->cw(edge))->info().is_in_interior() && !face->neighbor(face->cw(edge))->info().been_reconstructed()) {
		face->neighbor(face->cw(edge))->info().been_reconstructed(true);
    std::list<Triangulation::Vertex_handle> v1;
    get_boundary(face->neighbor(face->cw(edge)), face->neighbor(face->cw(edge))->index(face), v1);
		out_vertices.splice(out_vertices.end(), v1);
	}
	
	// Add central vertex
  out_vertices.push_back(face->vertex(edge));
	
	// Check counterclockwise edge
  if (face->neighbor(face->ccw(edge))->info().is_in_interior() && !face->neighbor(face->ccw(edge))->info().been_reconstructed()) {
		face->neighbor(face->ccw(edge))->info().been_reconstructed(true);
		std::list<Triangulation::Vertex_handle> v2;
    get_boundary(face->neighbor(face->ccw(edge)), face->neighbor(face->ccw(edge))->index(face), v2);
		out_vertices.splice(out_vertices.end(), v2);
	}
}