// prepair
//
// Copyright © 2009-2022,
// Ken Arroyo Ohori    k.ohori@tudelft.nl
// Hugo Ledoux         h.ledoux@tudelft.nl
// Martijn Meijers     b.m.meijers@tudelft.nl
// All rights reserved.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef Polygon_repair_h
#define Polygon_repair_h

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include "Enhanced_constrained_triangulation_2.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Exact_predicates_tag Tag;
struct VertexInfo {
  bool has_point;
  Kernel::Point_3 point;
  VertexInfo() {
    has_point = false;
    point = CGAL::ORIGIN;
  }
};
struct FaceInfo {
  bool processed;
  bool interior;
  FaceInfo() {
    processed = false;
    interior = false;
  }
};
typedef CGAL::Triangulation_vertex_base_with_info_2<VertexInfo, Kernel> VertexBase;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel> FaceBase;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo, Kernel, FaceBase> FaceBaseWithInfo;
typedef CGAL::Triangulation_data_structure_2<VertexBase, FaceBaseWithInfo> TriangulationDataStructure;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TriangulationDataStructure, Tag> ConstrainedDelaunayTriangulation;
typedef Enhanced_constrained_triangulation_2<ConstrainedDelaunayTriangulation> Triangulation;

struct Polygon_repair {
  OGRGeometry *geometry;
  Triangulation::Face_handle walk_start_location;
  Triangulation triangulation;
  std::list<Kernel::Point_3> points_in_polygon;
  Kernel::Plane_3 best_plane;
  
  void insert_points_in_list(OGRGeometry *g) {
    switch (g->getGeometryType()) {
        
      case wkbLineString25D: {
        OGRLinearRing *ring = static_cast<OGRLinearRing *>(g);
        ring->closeRings();
        for (int current_point = 1; current_point < ring->getNumPoints(); ++current_point) {
          points_in_polygon.push_back(Kernel::Point_3(ring->getX(current_point), ring->getY(current_point), ring->getZ(current_point)));
        } break;
      }
        
      case wkbPolygon25D: {
        OGRPolygon *polygon = static_cast<OGRPolygon *>(g);
        insert_points_in_list(polygon->getExteriorRing());
        for (int current_ring = 0; current_ring < polygon->getNumInteriorRings(); ++current_ring) insert_points_in_list(polygon->getInteriorRing(current_ring));
        break;
      }
        
      case wkbMultiPolygon25D: {
        OGRMultiPolygon *multipolygon = static_cast<OGRMultiPolygon *>(g);
        for (int current_polygon = 0; current_polygon < multipolygon->getNumGeometries(); ++current_polygon) {
          insert_points_in_list(multipolygon->getGeometryRef(current_polygon));
        } break;
      }
        
      default:
        std::cerr << "Error: input type << " << g->getGeometryType() << " << not supported" << std::endl;
        break;
        
    }
  }
  
  void compute_plane() {
    insert_points_in_list(geometry);
    linear_least_squares_fitting_3(points_in_polygon.begin(), points_in_polygon.end(), best_plane, CGAL::Dimension_tag<0>());
  }
  
  void insert_constraints_in_triangulation(OGRGeometry *g) {
    switch (g->getGeometryType()) {
      case wkbLineString: {
        OGRLinearRing *ring = static_cast<OGRLinearRing *>(g);
        ring->closeRings();
        Triangulation::Vertex_handle va, vb;
        vb = triangulation.insert(Kernel::Point_2(ring->getX(0), ring->getY(0)), walk_start_location);
        walk_start_location = triangulation.incident_faces(vb);
        for (int current_point = 1; current_point < ring->getNumPoints(); ++current_point) {
          va = vb;
          vb = triangulation.insert(Kernel::Point_2(ring->getX(current_point), ring->getY(current_point)), walk_start_location);
          if (va != vb) triangulation.odd_even_insert_constraint(va, vb);
          walk_start_location = triangulation.incident_faces(vb);
        } break;
      }
        
      case wkbLineString25D: {
        OGRLinearRing *ring = static_cast<OGRLinearRing *>(g);
        ring->closeRings();
        Triangulation::Vertex_handle va, vb;
        Kernel::Point_3 new_point(ring->getX(0), ring->getY(0), ring->getZ(0));
        vb = triangulation.insert(best_plane.to_2d(new_point), walk_start_location);
        vb->info().has_point = true;
        vb->info().point = new_point;
        walk_start_location = triangulation.incident_faces(vb);
        for (int current_point = 1; current_point < ring->getNumPoints(); ++current_point) {
          va = vb;
          Kernel::Point_3 new_point(ring->getX(current_point), ring->getY(current_point), ring->getZ(current_point));
          vb = triangulation.insert(best_plane.to_2d(new_point), walk_start_location);
          vb->info().has_point = true;
          vb->info().point = new_point;
          if (va != vb) triangulation.odd_even_insert_constraint(va, vb);
          walk_start_location = triangulation.incident_faces(vb);
        } break;
      }
        
      case wkbPolygon:
      case wkbPolygon25D: {
        OGRPolygon *polygon = static_cast<OGRPolygon *>(g);
        insert_constraints_in_triangulation(polygon->getExteriorRing());
        for (int current_ring = 0; current_ring < polygon->getNumInteriorRings(); ++current_ring) insert_constraints_in_triangulation(polygon->getInteriorRing(current_ring));
        break;
      }
        
      case wkbMultiPolygon:
      case wkbMultiPolygon25D: {
        OGRMultiPolygon *multipolygon = static_cast<OGRMultiPolygon *>(g);
        for (int current_polygon = 0; current_polygon < multipolygon->getNumGeometries(); ++current_polygon) {
          insert_constraints_in_triangulation(multipolygon->getGeometryRef(current_polygon));
        } break;
      }
        
      default:
        std::cerr << "Error: input type << " << g->getGeometryType() << " << not supported" << std::endl;
        break;
    }
  }
  
  void label_triangles() {
    std::list<Triangulation::Face_handle> to_check;
    triangulation.infinite_face()->info().processed = true;
    CGAL_assertion(triangulation.infinite_face()->info().processed == true);
    CGAL_assertion(triangulation.infinite_face()->info().interior == false);
    to_check.push_back(triangulation.infinite_face());
    while (!to_check.empty()) {
      CGAL_assertion(to_check.front()->info().processed == true);
      for (int neighbour = 0; neighbour < 3; ++neighbour) {
        if (to_check.front()->neighbor(neighbour)->info().processed == true) {
          // Note: validation code.
//          if (triangulation.is_constrained(Triangulation::Edge(to_check.front(), neighbour))) CGAL_assertion(to_check.front()->neighbor(neighbour)->info().interior != to_check.front()->info().interior);
//          else CGAL_assertion(to_check.front()->neighbor(neighbour)->info().interior == to_check.front()->info().interior);
        } else {
          to_check.front()->neighbor(neighbour)->info().processed = true;
          CGAL_assertion(to_check.front()->neighbor(neighbour)->info().processed == true);
          if (triangulation.is_constrained(Triangulation::Edge(to_check.front(), neighbour))) {
            to_check.front()->neighbor(neighbour)->info().interior = !to_check.front()->info().interior;
            to_check.push_back(to_check.front()->neighbor(neighbour));
          } else {
            to_check.front()->neighbor(neighbour)->info().interior = to_check.front()->info().interior;
            to_check.push_back(to_check.front()->neighbor(neighbour));
          }
        }
      } to_check.pop_front();
    }
  }
  
  void get_boundary(Triangulation::Face_handle face, int edge, std::list<Triangulation::Vertex_handle> &out_vertices) {
    // Check clockwise edge
    if (face->neighbor(face->cw(edge))->info().interior && !face->neighbor(face->cw(edge))->info().processed) {
      face->neighbor(face->cw(edge))->info().processed = true;
      std::list<Triangulation::Vertex_handle> v1;
      get_boundary(face->neighbor(face->cw(edge)), face->neighbor(face->cw(edge))->index(face), v1);
      out_vertices.splice(out_vertices.end(), v1);
    }
    
    // Add central vertex
    out_vertices.push_back(face->vertex(edge));
    
    // Check counterclockwise edge
    if (face->neighbor(face->ccw(edge))->info().interior && !face->neighbor(face->ccw(edge))->info().processed) {
      face->neighbor(face->ccw(edge))->info().processed = true;
      std::list<Triangulation::Vertex_handle> v2;
      get_boundary(face->neighbor(face->ccw(edge)), face->neighbor(face->ccw(edge))->index(face), v2);
      out_vertices.splice(out_vertices.end(), v2);
    }
  }
  
  void reconstruct() {
    if (triangulation.number_of_faces() < 1) {
      geometry = new OGRPolygon();
    }
    
    for (Triangulation::All_faces_iterator current_face = triangulation.all_faces_begin(); current_face != triangulation.all_faces_end(); ++current_face) {
      current_face->info().processed = false;
    } for (Triangulation::Finite_vertices_iterator current_vertex = triangulation.finite_vertices_begin(); current_vertex != triangulation.finite_vertices_end(); ++current_vertex) {
      if (!current_vertex->info().has_point) {
        current_vertex->info().point = best_plane.to_3d(current_vertex->point());
        current_vertex->info().has_point = true;
      }
    }
    
    // Reconstruct
    OGRMultiPolygon *out_geometries = new OGRMultiPolygon();
    for (Triangulation::Finite_faces_iterator seeding_face = triangulation.finite_faces_begin(); seeding_face != triangulation.finite_faces_end(); ++seeding_face) {
      
      if (!seeding_face->info().interior || seeding_face->info().processed) continue;
      seeding_face->info().processed = true;
      if (!seeding_face->info().processed) {
        std::cout << "Error: Face should be marked as reconstructed!" << std::endl;
      }
      
      // Get boundary
      std::list<Triangulation::Vertex_handle> vertices = std::list<Triangulation::Vertex_handle>();
      if (seeding_face->neighbor(2)->info().interior && !seeding_face->neighbor(2)->info().processed) {
        seeding_face->neighbor(2)->info().processed = true;
        std::list<Triangulation::Vertex_handle> l2;
        get_boundary(seeding_face->neighbor(2), seeding_face->neighbor(2)->index(seeding_face), l2);
        vertices.splice(vertices.end(), l2);
      } vertices.push_back(seeding_face->vertex(0));
      if (seeding_face->neighbor(1)->info().interior && !seeding_face->neighbor(1)->info().processed) {
        seeding_face->neighbor(1)->info().processed = true;
        std::list<Triangulation::Vertex_handle> l1;
        get_boundary(seeding_face->neighbor(1), seeding_face->neighbor(1)->index(seeding_face), l1);
        vertices.splice(vertices.end(), l1);
      } vertices.push_back(seeding_face->vertex(2));
      if (seeding_face->neighbor(0)->info().interior && !seeding_face->neighbor(0)->info().processed) {
        seeding_face->neighbor(0)->info().processed = true;
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
      std::list<std::list<Triangulation::Vertex_handle>> rings;
      std::stack<std::list<Triangulation::Vertex_handle>> chains_stack;
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
      for (std::list<std::list<Triangulation::Vertex_handle>>::iterator current_ring = rings.begin(); current_ring != rings.end(); ++current_ring) {
        std::list<Triangulation::Vertex_handle>::iterator smallest_vertex = current_ring->begin();
        for (std::list<Triangulation::Vertex_handle>::iterator current_vertex = current_ring->begin(); current_vertex != current_ring->end(); ++current_vertex) {
          if (geometry->Is3D()) {
            if ((*current_vertex)->info().point < (*smallest_vertex)->info().point) smallest_vertex = current_vertex;
          } else {
            if ((*current_vertex)->point() < (*smallest_vertex)->point()) smallest_vertex = current_vertex;
          }
         
        } if (current_ring->back() != *smallest_vertex) {
          ++smallest_vertex;
          current_ring->splice(current_ring->begin(), *current_ring, smallest_vertex, current_ring->end());
        }
      }
      
      // Make rings
      if (rings.size() == 0) continue;
      std::list<OGRLinearRing *> rings_for_polygon;
      for (std::list<std::list<Triangulation::Vertex_handle>>::iterator current_ring = rings.begin(); current_ring != rings.end(); ++current_ring) {
        OGRLinearRing *new_ring = new OGRLinearRing();
        for (std::list<Triangulation::Vertex_handle>::reverse_iterator current_vertex = current_ring->rbegin(); current_vertex != current_ring->rend(); ++current_vertex) {
          if (geometry->Is3D()) new_ring->addPoint((*current_vertex)->info().point.x(),
                                                   (*current_vertex)->info().point.y(),
                                                   (*current_vertex)->info().point.z());
          else new_ring->addPoint((*current_vertex)->point().x(), (*current_vertex)->point().y());
        } if (geometry->Is3D()) new_ring->addPoint(current_ring->back()->info().point.x(),
                                                   current_ring->back()->info().point.y(),
                                                   current_ring->back()->info().point.z());
        else new_ring->addPoint(current_ring->back()->point().x(), current_ring->back()->point().y());
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
      geometry = new OGRPolygon();
    }
    
    else if (out_geometries->getNumGeometries() == 1) {
      OGRPolygon *new_polygon = static_cast<OGRPolygon *>(out_geometries->getGeometryRef(0)->clone());
      delete out_geometries;
      geometry = new_polygon;
    }
    
    else geometry = out_geometries;
  }
  
  void repair() {
    if (geometry->Is3D()) compute_plane();
    insert_constraints_in_triangulation(geometry);
    if (triangulation.number_of_faces() == 0) {
      geometry = NULL;
      return;
    } label_triangles();
    reconstruct();
  }
};

#endif /* Polygon_repair_h */
