/*
 Copyright (c) 2009-2013,
 Ken Arroyo Ohori    g.a.k.arroyoohori@tudelft.nl
 Hugo Ledoux         h.ledoux@tudelft.nl
 Martijn Meijers     b.m.meijers@tudelft.nl
 All rights reserved.
 
 // This file is part of prepair: you can redistribute it and/or modify
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

#include "PolygonRepair.h"

void PolygonRepair::insertOuterConstraints(Triangulation &triangulation, OGRGeometry* geometry) {
    Triangulation::Vertex_handle va, vb;
    Triangulation::Face_handle faceOfEdge;
    int indexOfEdge;
    
    triangulation.clear();
    
    switch (geometry->getGeometryType()) {
        case wkbPolygon: {
            OGRPolygon *polygon = (OGRPolygon *)geometry;
            
            for (int currentPoint = 0; currentPoint < polygon->getExteriorRing()->getNumPoints(); ++currentPoint) {
                va = triangulation.insert(Point(polygon->getExteriorRing()->getX(currentPoint),
                                                polygon->getExteriorRing()->getY(currentPoint)));
                vb = triangulation.insert(Point(polygon->getExteriorRing()->getX((currentPoint+1)%polygon->getExteriorRing()->getNumPoints()),
                                                polygon->getExteriorRing()->getY((currentPoint+1)%polygon->getExteriorRing()->getNumPoints())));
                if (triangulation.is_edge(va, vb, faceOfEdge, indexOfEdge)) {
                    if (triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(faceOfEdge, indexOfEdge))) {
                        std::cerr << "Removing duplicate constraint <" << va->point() << ", " << vb->point() << "> in outer ring..." << std::endl;
                        triangulation.insert_constraint(va, vb);
                        triangulation.remove_constraint(va, vb);
                    } else triangulation.insert_constraint(va, vb);
                } else triangulation.insert_constraint(va, vb);
            } break;
        }
            
        case wkbMultiPolygon: {
            OGRMultiPolygon *multipolygon = static_cast<OGRMultiPolygon *>(geometry);
            for (int currentPolygon = 0; currentPolygon < multipolygon->getNumGeometries(); ++currentPolygon) {
                OGRPolygon *polygon = static_cast<OGRPolygon *>(multipolygon->getGeometryRef(currentPolygon));
                
                // Outer
                for (int currentPoint = 0; currentPoint < polygon->getExteriorRing()->getNumPoints(); ++currentPoint) {
                    va = triangulation.insert(Point(polygon->getExteriorRing()->getX(currentPoint),
                                                    polygon->getExteriorRing()->getY(currentPoint)));
                    vb = triangulation.insert(Point(polygon->getExteriorRing()->getX((currentPoint+1)%polygon->getExteriorRing()->getNumPoints()),
                                                    polygon->getExteriorRing()->getY((currentPoint+1)%polygon->getExteriorRing()->getNumPoints())));
                    if (triangulation.is_edge(va, vb, faceOfEdge, indexOfEdge)) {
                        if (triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(faceOfEdge, indexOfEdge))) {
                            std::cerr << "Removing duplicate constraint <" << va->point() << ", " << vb->point() << "> in outer ring..." << std::endl;
                            triangulation.insert_constraint(va, vb);
                            triangulation.remove_constraint(va, vb);
                        } else triangulation.insert_constraint(va, vb);
                    } else triangulation.insert_constraint(va, vb);
                }
                
                // Inner
                for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing) {
                    for (int currentPoint = 0; currentPoint < polygon->getInteriorRing(currentRing)->getNumPoints(); ++currentPoint) {
                        va = triangulation.insert(Point(polygon->getInteriorRing(currentRing)->getX(currentPoint),
                                                        polygon->getInteriorRing(currentRing)->getY(currentPoint)));
                        vb = triangulation.insert(Point(polygon->getInteriorRing(currentRing)->getX((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints()),
                                                        polygon->getInteriorRing(currentRing)->getY((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints())));
                        if (triangulation.is_edge(va, vb, faceOfEdge, indexOfEdge)) {
                            if (triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(faceOfEdge, indexOfEdge))) {
                                std::cerr << "Removing duplicate constraint <" << va->point() << ", " << vb->point() << "> in outer ring..." << std::endl;
                                triangulation.insert_constraint(va, vb);
                                triangulation.remove_constraint(va, vb);
                            } else triangulation.insert_constraint(va, vb);
                        } else triangulation.insert_constraint(va, vb);
                    }
                }
            } break;
        }
        
        default:
            std::cerr << "Error: Cannot understand input. Only polygons are supported." << std::endl;
            break;
    }
    
    // Remove partially even-overlapping constraints
    for (Triangulation::Subconstraint_iterator currentEdge = triangulation.subconstraints_begin();
         currentEdge != triangulation.subconstraints_end();
         ++currentEdge) {
        if (triangulation.number_of_enclosing_constraints(currentEdge->first.first, currentEdge->first.second) % 2 == 0) {
            std::cerr << "Removing even-overlapped subconstraint <" << currentEdge->first.first->point() << ", " << currentEdge->first.second->point() << ">..." << std::endl;
            triangulation.remove_constraint(currentEdge->first.first, currentEdge->first.second);
        }
    }
}

void PolygonRepair::insertConstraints(Triangulation &triangulation, OGRGeometry* geometry) {
    Triangulation::Vertex_handle va, vb;
    Triangulation::Face_handle faceOfEdge;
    int indexOfEdge;
    
    triangulation.clear();
    
    switch (geometry->getGeometryType()) {
            
        case wkbPolygon: {
            OGRPolygon *polygon = static_cast<OGRPolygon *>(geometry);
            
            // Outer
            for (int currentPoint = 0; currentPoint < polygon->getExteriorRing()->getNumPoints(); ++currentPoint) {
                va = triangulation.insert(Point(polygon->getExteriorRing()->getX(currentPoint),
                                                polygon->getExteriorRing()->getY(currentPoint)));
                vb = triangulation.insert(Point(polygon->getExteriorRing()->getX((currentPoint+1)%polygon->getExteriorRing()->getNumPoints()),
                                                polygon->getExteriorRing()->getY((currentPoint+1)%polygon->getExteriorRing()->getNumPoints())));
                if (triangulation.is_edge(va, vb, faceOfEdge, indexOfEdge)) {
                    if (triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(faceOfEdge, indexOfEdge))) {
                        std::cerr << "Removing duplicate constraint <" << va->point() << ", " << vb->point() << "> in outer ring..." << std::endl;
                        triangulation.insert_constraint(va, vb);
                        triangulation.remove_constraint(va, vb);
                    } else triangulation.insert_constraint(va, vb);
                } else triangulation.insert_constraint(va, vb);
            }
            
            // Inner
            for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing) {
                for (int currentPoint = 0; currentPoint < polygon->getInteriorRing(currentRing)->getNumPoints(); ++currentPoint) {
                    va = triangulation.insert(Point(polygon->getInteriorRing(currentRing)->getX(currentPoint),
                                                    polygon->getInteriorRing(currentRing)->getY(currentPoint)));
                    vb = triangulation.insert(Point(polygon->getInteriorRing(currentRing)->getX((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints()),
                                                    polygon->getInteriorRing(currentRing)->getY((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints())));
                    if (triangulation.is_edge(va, vb, faceOfEdge, indexOfEdge)) {
                        if (triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(faceOfEdge, indexOfEdge))) {
                            std::cerr << "Removing duplicate constraint <" << va->point() << ", " << vb->point() << "> in outer ring..." << std::endl;
                            triangulation.insert_constraint(va, vb);
                            triangulation.remove_constraint(va, vb);
                        } else triangulation.insert_constraint(va, vb);
                    } else triangulation.insert_constraint(va, vb);
                }
            } break;
            
        }
        
        case wkbMultiPolygon: {
            OGRMultiPolygon *multipolygon = static_cast<OGRMultiPolygon *>(geometry);
            for (int currentPolygon = 0; currentPolygon < multipolygon->getNumGeometries(); ++currentPolygon) {
                OGRPolygon *polygon = static_cast<OGRPolygon *>(multipolygon->getGeometryRef(currentPolygon));
                
                // Outer
                for (int currentPoint = 0; currentPoint < polygon->getExteriorRing()->getNumPoints(); ++currentPoint) {
                    va = triangulation.insert(Point(polygon->getExteriorRing()->getX(currentPoint),
                                                    polygon->getExteriorRing()->getY(currentPoint)));
                    vb = triangulation.insert(Point(polygon->getExteriorRing()->getX((currentPoint+1)%polygon->getExteriorRing()->getNumPoints()),
                                                    polygon->getExteriorRing()->getY((currentPoint+1)%polygon->getExteriorRing()->getNumPoints())));
                    if (triangulation.is_edge(va, vb, faceOfEdge, indexOfEdge)) {
                        if (triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(faceOfEdge, indexOfEdge))) {
                            std::cerr << "Removing duplicate constraint <" << va->point() << ", " << vb->point() << "> in outer ring..." << std::endl;
                            triangulation.insert_constraint(va, vb);
                            triangulation.remove_constraint(va, vb);
                        } else triangulation.insert_constraint(va, vb);
                    } else triangulation.insert_constraint(va, vb);
                }
                
                // Inner
                for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing) {
                    for (int currentPoint = 0; currentPoint < polygon->getInteriorRing(currentRing)->getNumPoints(); ++currentPoint) {
                        va = triangulation.insert(Point(polygon->getInteriorRing(currentRing)->getX(currentPoint),
                                                        polygon->getInteriorRing(currentRing)->getY(currentPoint)));
                        vb = triangulation.insert(Point(polygon->getInteriorRing(currentRing)->getX((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints()),
                                                        polygon->getInteriorRing(currentRing)->getY((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints())));
                        if (triangulation.is_edge(va, vb, faceOfEdge, indexOfEdge)) {
                            if (triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(faceOfEdge, indexOfEdge))) {
                                std::cerr << "Removing duplicate constraint <" << va->point() << ", " << vb->point() << "> in outer ring..." << std::endl;
                                triangulation.insert_constraint(va, vb);
                                triangulation.remove_constraint(va, vb);
                            } else triangulation.insert_constraint(va, vb);
                        } else triangulation.insert_constraint(va, vb);
                    }
                }
            } break;
        }
        
        default:
            std::cerr << "Error: Cannot understand input. Only polygons are supported." << std::endl;
            break;
    }
    
    // Remove partially even-overlapping constraints
    for (Triangulation::Subconstraint_iterator currentEdge = triangulation.subconstraints_begin();
         currentEdge != triangulation.subconstraints_end();
         ++currentEdge) {
        if (triangulation.number_of_enclosing_constraints(currentEdge->first.first, currentEdge->first.second) % 2 == 0) {
            std::cerr << "Removing even-overlapped subconstraint <" << currentEdge->first.first->point() << ", " << currentEdge->first.second->point() << ">..." << std::endl;
            triangulation.remove_constraint(currentEdge->first.first, currentEdge->first.second);
        }
    }
}

void PolygonRepair::tagOddEven(Triangulation &triangulation) {
	
    // Clean tags
    for (Triangulation::Face_handle currentFace = triangulation.all_faces_begin(); currentFace != triangulation.all_faces_end(); ++currentFace)
        currentFace->info().clear();
    
    // Initialise tagging
    std::stack<Triangulation::Face_handle> interiorStack, exteriorStack;
    exteriorStack.push(triangulation.infinite_face());
    std::stack<Triangulation::Face_handle> *currentStack = &exteriorStack;
    std::stack<Triangulation::Face_handle> *dualStack = &interiorStack;
    bool labellingInterior = false;
    
    
    // Until we finish
    while (!interiorStack.empty() || !exteriorStack.empty()) {
        
        // Give preference to whatever we're already doing
        while (!currentStack->empty()) {
            Triangulation::Face_handle currentFace = currentStack->top();
			currentStack->pop();
            if (currentFace->info().hasBeenProcessed()) continue;
			currentFace->info().isInInterior(labellingInterior);
            for (int currentEdge = 0; currentEdge < 3; ++currentEdge) {
                if (!currentFace->neighbor(currentEdge)->info().hasBeenProcessed()) {
                    if (currentFace->is_constrained(currentEdge))
                        dualStack->push(currentFace->neighbor(currentEdge));
                    else
                        currentStack->push(currentFace->neighbor(currentEdge));
                }
            }
        }
        
        // Flip
        if (!labellingInterior) {
            currentStack = &interiorStack;
            dualStack = &exteriorStack;
        } else {
            currentStack = &exteriorStack;
            dualStack = &interiorStack;
        } labellingInterior = !labellingInterior;
	}
}

void PolygonRepair::tagPointSet(Triangulation &triangulation, OGRGeometry* geometry) {
    
    // Clean tags
    for (Triangulation::Face_handle currentFace = triangulation.all_faces_begin(); currentFace != triangulation.all_faces_end(); ++currentFace)
        currentFace->info().clear();
    
    // Tag everything connected to the infinite face plus the triangles adjancent to those
    switch (geometry->getGeometryType()) {
            
        case wkbPolygon: {
            OGRPolygon *polygon = (OGRPolygon *)geometry;
            for (int currentPoint = 0; currentPoint < polygon->getExteriorRing()->getNumPoints(); ++currentPoint) {
                
            }
            break;
        }
            
        default:
            std::cerr << "Error: Cannot understand input. Only polygons are supported." << std::endl;
            break;
            
    }
}

void PolygonRepair::printEdges(Triangulation &triangulation) {
    std::cout << "Edges:" << std::endl;
    for (Triangulation::Finite_edges_iterator currentEdge = triangulation.finite_edges_begin();
         currentEdge != triangulation.finite_edges_end();
         ++currentEdge) {
        if (triangulation.is_constrained(*currentEdge))
            std::cout << "\t<" << currentEdge->first->vertex((currentEdge->second+1)%3)->point() << ", " << currentEdge->first->vertex((currentEdge->second+2)%3)->point() << ">" << std::endl;
    }
}