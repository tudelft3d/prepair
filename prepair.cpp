/*
 Copyright (c) 2009-2013,
 Ken Arroyo Ohori    g.a.k.arroyoohori@tudelft.nl
 Hugo Ledoux         h.ledoux@tudelft.nl
 Martijn Meijers     b.m.meijers@tudelft.nl
 All rights reserved.
 
 This file is part of pprepair: you can redistribute it and/or modify
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
// if the code crashes, try to compile with EXACT_CONSTRUCTIONS so that
// robust arithmetic is used
#define EXACT_CONSTRUCTIONS

// STL
#include <iostream>
#include <stack>
#include <set>
#include <fstream>

// OGR
#include <gdal/ogrsf_frmts.h>

// CGAL
#ifdef EXACT_CONSTRUCTIONS
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#else
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#endif
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

// Kernel
#ifdef EXACT_CONSTRUCTIONS
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
#else
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
#endif

typedef CGAL::Triangulation_vertex_base_2<K> VB;
typedef CGAL::Constrained_triangulation_face_base_2<K> FB;
typedef CGAL::Triangulation_face_base_with_info_2<void *, K, FB> FBWI;
typedef CGAL::Triangulation_data_structure_2<VB, FBWI> TDS;
typedef CGAL::Exact_predicates_tag PT;
typedef CGAL::Exact_intersections_tag IT;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, PT> CDT;

typedef CGAL::Constrained_triangulation_plus_2<CDT> Triangulation;
typedef Triangulation::Point Point;


OGRMultiPolygon* repair(OGRGeometry* geometry);
void tag(Triangulation &triangulation, void *interiorHandle, void *exteriorHandle);
std::list<Triangulation::Vertex_handle> *getBoundary(Triangulation::Face_handle face, int edge);




int main (int argc, const char * argv[]) {
  
  if (argc < 2 || argc > 3 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
    std::cout << "=== prepair Help ===\n" << std::endl;
    std::cout << "Usage:   prepair 'POLYGON(...)'" << std::endl;
    std::cout << "OR" << std::endl;
    std::cout << "Usage:   prepair -f infile.txt (infile.txt must contain one WKT on the 1st line)" << std::endl;
    return 0;
  }
  
  // Read input
  unsigned int bufferSize = 10000000;
  char *inputWKT = (char *)malloc(bufferSize*sizeof(char *));
  
  for (int argNum = 1; argNum < argc; ++argNum) {
    if (strcmp(argv[argNum], "-f") == 0) {
      
      if (argNum + 1 <= argc - 1 && argv[argNum+1][0] != '-') {
        std::ifstream infile(argv[argNum+1], std::ifstream::in);
        infile.getline(inputWKT, bufferSize);
        ++argNum;
      } else {
        std::cerr << "Error: Missing input file name." << std::endl;
        return 1;
      }
    }
    else 
      strcpy(inputWKT, argv[argNum]);
  }
  
  std::cout << "Processing: " << inputWKT << std::endl;
  
  OGRGeometry *geometry;
  OGRGeometryFactory::createFromWkt(&inputWKT, NULL, &geometry);
  if (geometry == NULL) {
    std::cout << "Error: WKT is not valid" << std::endl;
    return 1;
  }
  
  OGRMultiPolygon* outputPolygons = repair(geometry);
  
  if (outputPolygons == NULL) {
    std::cout << "Impossible to repair the polygon: input points are collinear (no area given)." << std::endl;
    return 0;
  }
  else {
    char *outputWKT;
    outputPolygons->exportToWkt(&outputWKT);
    std::cout << std::endl << "Repaired polygon:" << std::endl << outputWKT << std::endl;
    return 0;
  }
}



void tag(Triangulation &triangulation, void *interiorHandle, void *exteriorHandle) {
	
    // Clean tags
    for (Triangulation::Face_handle currentFace = triangulation.all_faces_begin(); currentFace != triangulation.all_faces_end(); ++currentFace)
        currentFace->info() = NULL;
    
    // Initialise tagging
    std::stack<Triangulation::Face_handle> interiorStack, exteriorStack;
    exteriorStack.push(triangulation.infinite_face());
    std::stack<Triangulation::Face_handle> *currentStack = &exteriorStack;
    std::stack<Triangulation::Face_handle> *dualStack = &interiorStack;
    void *currentHandle = exteriorHandle;
    void *dualHandle = interiorHandle;
    
    // Until we finish
    while (!interiorStack.empty() || !exteriorStack.empty()) {
        
        // Give preference to whatever we're already doing
        while (!currentStack->empty()) {
            Triangulation::Face_handle currentFace = currentStack->top();
			currentStack->pop();
            if (currentFace->info() != NULL) continue;
			currentFace->info() = currentHandle;
            for (int currentEdge = 0; currentEdge < 3; ++currentEdge) {
                if (currentFace->neighbor(currentEdge)->info() == NULL)
                    if (currentFace->is_constrained(currentEdge)) dualStack->push(currentFace->neighbor(currentEdge));
                else currentStack->push(currentFace->neighbor(currentEdge));
            }
        }
			
        // Flip
        if (currentHandle == exteriorHandle) {
            currentHandle = interiorHandle;
            dualHandle = exteriorHandle;
            currentStack = &interiorStack;
            dualStack = &exteriorStack;
        } else {
            currentHandle = exteriorHandle;
            dualHandle = interiorHandle;
            currentStack = &exteriorStack;
            dualStack = &interiorStack;
        }
	}
}




std::list<Triangulation::Vertex_handle> *getBoundary(Triangulation::Face_handle face, int edge) {
    
    std::list<Triangulation::Vertex_handle> *vertices = new std::list<Triangulation::Vertex_handle>();
    
    // Check clockwise edge
    if (!face->is_constrained(face->cw(edge)) && face->neighbor(face->cw(edge))->info() != NULL) {
		face->neighbor(face->cw(edge))->info() = NULL;
		std::list<Triangulation::Vertex_handle> *v1 = getBoundary(face->neighbor(face->cw(edge)), face->neighbor(face->cw(edge))->index(face));
		vertices->splice(vertices->end(), *v1);
		delete v1;
	}
	
	// Add central vertex
	vertices->push_back(face->vertex(edge));
	
	// Check counterclockwise edge
    if (!face->is_constrained(face->ccw(edge)) && face->neighbor(face->ccw(edge))->info() != NULL) {
		face->neighbor(face->ccw(edge))->info() = NULL;
		std::list<Triangulation::Vertex_handle> *v2 = getBoundary(face->neighbor(face->ccw(edge)), face->neighbor(face->ccw(edge))->index(face));
		vertices->splice(vertices->end(), *v2);
		delete v2;
	}
	
    return vertices;
}




OGRMultiPolygon* repair(OGRGeometry* geometry) {
  
  // Triangulation
  Triangulation triangulation;
  switch (geometry->getGeometryType()) {
      
    case wkbPolygon: {
      
      OGRPolygon *polygon = (OGRPolygon *)geometry;
      for (int currentPoint = 0; currentPoint < polygon->getExteriorRing()->getNumPoints(); ++currentPoint) {
        triangulation.insert_constraint(Point(polygon->getExteriorRing()->getX(currentPoint), 
                                              polygon->getExteriorRing()->getY(currentPoint)),
                                        Point(polygon->getExteriorRing()->getX((currentPoint+1)%polygon->getExteriorRing()->getNumPoints()), 
                                              polygon->getExteriorRing()->getY((currentPoint+1)%polygon->getExteriorRing()->getNumPoints())));
      } for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing) {
        for (int currentPoint = 0; currentPoint < polygon->getInteriorRing(currentRing)->getNumPoints(); ++currentPoint)
          triangulation.insert_constraint(Point(polygon->getInteriorRing(currentRing)->getX(currentPoint), 
                                                polygon->getInteriorRing(currentRing)->getY(currentPoint)),
                                          Point(polygon->getInteriorRing(currentRing)->getX((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints()), 
                                                polygon->getInteriorRing(currentRing)->getY((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints())));
      } break;
      
    } default:
      std::cout << "Error: Cannot understand input. Only polygons are supported." << std::endl;
      break;
  } 
  
//  std::cout << "Triangulation: " << triangulation.number_of_faces() << " faces, " << triangulation.number_of_vertices() << " vertices." << std::endl;
  if (triangulation.number_of_faces() < 1) {
    return NULL;
  }
  
  // Tag
  void *interior = malloc(sizeof(void *));
  void *exterior = malloc(sizeof(void *));
  tag(triangulation, interior, exterior);
  
  // Reconstruct
  //    OGRMultiPolygon outputPolygons;
  OGRMultiPolygon* outputPolygons = new OGRMultiPolygon();
  for (Triangulation::Finite_faces_iterator seedingFace = triangulation.finite_faces_begin(); seedingFace != triangulation.finite_faces_end(); ++seedingFace) {
    
    if (seedingFace->info() != interior) continue;
    
    // Get boundary
    std::list<Triangulation::Vertex_handle> *vertices = new std::list<Triangulation::Vertex_handle>();
    seedingFace->info() = NULL;
    if (seedingFace->neighbor(2)->info() == interior) {
      seedingFace->neighbor(2)->info() = NULL;
      std::list<Triangulation::Vertex_handle> *l2 = getBoundary(seedingFace->neighbor(2), seedingFace->neighbor(2)->index(seedingFace));
      vertices->splice(vertices->end(), *l2);
      delete l2;
    } vertices->push_back(seedingFace->vertex(0));
    if (seedingFace->neighbor(1)->info() == interior) {
      seedingFace->neighbor(1)->info() = NULL;
      std::list<Triangulation::Vertex_handle> *l1 = getBoundary(seedingFace->neighbor(1), seedingFace->neighbor(1)->index(seedingFace));
      vertices->splice(vertices->end(), *l1);
      delete l1;
    } vertices->push_back(seedingFace->vertex(2));
    if (seedingFace->neighbor(0)->info() == interior) {
      seedingFace->neighbor(0)->info() = NULL;
      std::list<Triangulation::Vertex_handle> *l0 = getBoundary(seedingFace->neighbor(0), seedingFace->neighbor(0)->index(seedingFace));
      vertices->splice(vertices->end(), *l0);
      delete l0;
    } vertices->push_back(seedingFace->vertex(1));
    
    // Find cutting vertices
    std::set<Triangulation::Vertex_handle> visitedVertices;
    std::set<Triangulation::Vertex_handle> repeatedVertices;
    for (std::list<Triangulation::Vertex_handle>::iterator currentVertex = vertices->begin(); currentVertex != vertices->end(); ++currentVertex) {
      if (!visitedVertices.insert(*currentVertex).second) repeatedVertices.insert(*currentVertex);
    } visitedVertices.clear();
    
    // Cut and join rings in the correct order
    std::list<std::list<Triangulation::Vertex_handle> *> rings;
    std::stack<std::list<Triangulation::Vertex_handle> *> chainsStack;
    std::map<Triangulation::Vertex_handle, std::list<Triangulation::Vertex_handle> *> vertexChainMap;
    std::list<Triangulation::Vertex_handle> *newChain = new std::list<Triangulation::Vertex_handle>();
    for (std::list<Triangulation::Vertex_handle>::iterator currentVertex = vertices->begin(); currentVertex != vertices->end(); ++currentVertex) {
      
      // New chain
      if (repeatedVertices.count(*currentVertex) > 0) {
        // Closed by itself
        if (newChain->front() == *currentVertex) {
          // Degenerate (insufficient vertices to be valid)
          if (newChain->size() < 3) delete newChain;
          else {
            std::list<Triangulation::Vertex_handle>::iterator secondElement = newChain->begin();
            ++secondElement;
            // Degenerate (zero area)
            if (newChain->back() == *secondElement) delete newChain;
            // Valid
            else rings.push_back(newChain);
          }
        }
        // Open by itself
        else {
          // Closed with others in stack
          if (vertexChainMap.count(*currentVertex)) {
            
            while (chainsStack.top() != vertexChainMap[*currentVertex]) {
              newChain->splice(newChain->begin(), *chainsStack.top());
              chainsStack.pop();
            } newChain->splice(newChain->begin(), *chainsStack.top());
            chainsStack.pop();
            vertexChainMap.erase(*currentVertex);
            // Degenerate (insufficient vertices to be valid)
            if (newChain->size() < 3) delete newChain;
            else {
              std::list<Triangulation::Vertex_handle>::iterator secondElement = newChain->begin();
              ++secondElement;
              // Degenerate (zero area)
              if (newChain->back() == *secondElement) delete newChain;
              // Valid
              else rings.push_back(newChain);
            }
          }
          // Open
          else {
            // Not first chain
            if (repeatedVertices.count(newChain->front()) > 0) vertexChainMap[newChain->front()] = newChain;
            chainsStack.push(newChain);
          }
        } newChain = new std::list<Triangulation::Vertex_handle>();
      } newChain->push_back(*currentVertex);
    }
    // Final ring
    while (chainsStack.size() > 0) {
      newChain->splice(newChain->begin(), *chainsStack.top());
      chainsStack.pop();
    }
    // Degenerate (insufficient vertices to be valid)
    if (newChain->size() < 3) delete newChain;
    else {
      std::list<Triangulation::Vertex_handle>::iterator secondElement = newChain->begin();
      ++secondElement;
      // Degenerate (zero area)
      if (newChain->back() == *secondElement) delete newChain;
      // Valid
      else rings.push_back(newChain);
    }
    // Make rings
    std::list<OGRLinearRing *> ringsForPolygon;
    for (std::list<std::list<Triangulation::Vertex_handle> *>::iterator currentRing = rings.begin(); currentRing != rings.end(); ++currentRing) {
      OGRLinearRing *newRing = new OGRLinearRing();
      for (std::list<Triangulation::Vertex_handle>::reverse_iterator currentVertex = (*currentRing)->rbegin(); currentVertex != (*currentRing)->rend(); ++currentVertex) {
        newRing->addPoint(CGAL::to_double((*currentVertex)->point().x()), CGAL::to_double((*currentVertex)->point().y()));
      } newRing->addPoint(CGAL::to_double((*currentRing)->back()->point().x()), CGAL::to_double((*currentRing)->back()->point().y()));
      ringsForPolygon.push_back(newRing);
    } OGRPolygon *newPolygon = new OGRPolygon();
    for (std::list<OGRLinearRing *>::iterator currentRing = ringsForPolygon.begin(); currentRing != ringsForPolygon.end(); ++currentRing) {
      if (!(*currentRing)->isClockwise()) {
        newPolygon->addRingDirectly(*currentRing);
        break;
      }
    } for (std::list<OGRLinearRing *>::iterator currentRing = ringsForPolygon.begin(); currentRing != ringsForPolygon.end(); ++currentRing)
      if ((*currentRing)->isClockwise()) newPolygon->addRingDirectly(*currentRing);
    outputPolygons->addGeometryDirectly(newPolygon);
  }
  return outputPolygons;
}

