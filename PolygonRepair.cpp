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

OGRMultiPolygon *PolygonRepair::repairOddEven(OGRGeometry *geometry) {
    Triangulation triangulation;
    insertConstraints(triangulation, geometry);
    tagOddEven(triangulation);
    OGRMultiPolygon *outPolygons = reconstruct(triangulation);
    return outPolygons;
}

OGRMultiPolygon *PolygonRepair::repairPointSet(OGRGeometry *geometry) {
    std::list<std::pair<bool, OGRMultiPolygon *> > repairedRings;   // bool indicates if outer/inner are flipped
    
    switch (geometry->getGeometryType()) {
            
        case wkbLinearRing: {
            repairedRings.push_back(std::pair<bool, OGRMultiPolygon *>(false, repairOddEven(geometry)));
            break;
        }
            
        case wkbPolygon: {
            OGRPolygon *polygon = static_cast<OGRPolygon *>(geometry);
            repairedRings.push_back(std::pair<bool, OGRMultiPolygon *>(false, repairOddEven(polygon->getExteriorRing())));
            for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing)
                repairedRings.push_back(std::pair<bool, OGRMultiPolygon *>(true, repairOddEven(polygon->getInteriorRing(currentRing))));
            break;
        }
            
        case wkbMultiPolygon: {
            OGRMultiPolygon *multipolygon = static_cast<OGRMultiPolygon *>(geometry);
            for (int currentPolygon = 0; currentPolygon < multipolygon->getNumGeometries(); ++currentPolygon) {
                OGRPolygon *polygon = static_cast<OGRPolygon *>(multipolygon->getGeometryRef(currentPolygon));
                repairedRings.push_back(std::pair<bool, OGRMultiPolygon *>(false, repairOddEven(polygon->getExteriorRing())));
                for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing)
                    repairedRings.push_back(std::pair<bool, OGRMultiPolygon *>(true, repairOddEven(polygon->getInteriorRing(currentRing))));
            } break;
        }
            
        default:
            std::cerr << "Error: Cannot understand input. Only polygons are supported." << std::endl;
            break;
    }
    
    Triangulation triangulation;
    for (std::list<std::pair<bool, OGRMultiPolygon *> >::iterator currentMultipolygon = repairedRings.begin();
         currentMultipolygon != repairedRings.end(); ++currentMultipolygon)
        insertConstraints(triangulation, currentMultipolygon->second);
    tagPointSet(triangulation, repairedRings);
    OGRMultiPolygon *outPolygons = reconstruct(triangulation);
    return outPolygons;
}

void PolygonRepair::insertConstraints(Triangulation &triangulation, OGRGeometry* geometry) {
    Triangulation::Vertex_handle va, vb;
    Triangulation::Face_handle faceOfEdge;
    int indexOfEdge;
    
    triangulation.clear();
    
    switch (geometry->getGeometryType()) {
            
        case wkbLinearRing: {
            OGRLinearRing *ring = static_cast<OGRLinearRing *>(geometry);
            
            for (int currentPoint = 0; currentPoint < ring->getNumPoints(); ++currentPoint) {
                va = triangulation.insert(Point(ring->getX(currentPoint),
                                                ring->getY(currentPoint)));
                vb = triangulation.insert(Point(ring->getX((currentPoint+1)%ring->getNumPoints()),
                                                ring->getY((currentPoint+1)%ring->getNumPoints())));
                std::cerr << "Inserting contraint <" << va->point() << ", " << vb->point() << "> in outer ring..." << std::endl;
                if (triangulation.is_edge(va, vb, faceOfEdge, indexOfEdge)) {
                    if (triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(faceOfEdge, indexOfEdge))) {
                        std::cerr << "Removing duplicate constraint <" << va->point() << ", " << vb->point() << "> in outer ring..." << std::endl;
                        triangulation.insert_constraint(va, vb);
                        triangulation.remove_constraint(va, vb);
                    } else triangulation.insert_constraint(va, vb);
                } else triangulation.insert_constraint(va, vb);
            }
        }
            
        case wkbPolygon: {
            OGRPolygon *polygon = static_cast<OGRPolygon *>(geometry);
            
            // Outer
            for (int currentPoint = 0; currentPoint < polygon->getExteriorRing()->getNumPoints(); ++currentPoint) {
                va = triangulation.insert(Point(polygon->getExteriorRing()->getX(currentPoint),
                                                polygon->getExteriorRing()->getY(currentPoint)));
                vb = triangulation.insert(Point(polygon->getExteriorRing()->getX((currentPoint+1)%polygon->getExteriorRing()->getNumPoints()),
                                                polygon->getExteriorRing()->getY((currentPoint+1)%polygon->getExteriorRing()->getNumPoints())));
                std::cerr << "Inserting contraint <" << va->point() << ", " << vb->point() << "> in outer ring..." << std::endl;
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
                    std::cerr << "Inserting contraint <" << va->point() << ", " << vb->point() << "> in inner ring..." << std::endl;
                    if (triangulation.is_edge(va, vb, faceOfEdge, indexOfEdge)) {
                        if (triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(faceOfEdge, indexOfEdge))) {
                            std::cerr << "Removing duplicate constraint <" << va->point() << ", " << vb->point() << "> in inner ring..." << std::endl;
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
                    std::cerr << "Inserting contraint <" << va->point() << ", " << vb->point() << "> in outer ring..." << std::endl;
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
                        std::cerr << "Inserting contraint <" << va->point() << ", " << vb->point() << "> in inner ring..." << std::endl;
                        if (triangulation.is_edge(va, vb, faceOfEdge, indexOfEdge)) {
                            if (triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(faceOfEdge, indexOfEdge))) {
                                std::cerr << "Removing duplicate constraint <" << va->point() << ", " << vb->point() << "> in inner ring..." << std::endl;
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
            if (currentFace->info().beenTagged()) continue;
			currentFace->info().isInInterior(labellingInterior);
            for (int currentEdge = 0; currentEdge < 3; ++currentEdge) {
                if (!currentFace->neighbor(currentEdge)->info().beenTagged()) {
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

void PolygonRepair::tagPointSet(Triangulation &triangulation, std::list<std::pair<bool, OGRMultiPolygon *> > &geometries) {
    std::stack<Triangulation::Face_handle> borderTriangles, taggingStack;
    Triangulation::Vertices_in_constraint_iterator currentVertex, nextVertex;
    
    for (std::list<std::pair<bool, OGRMultiPolygon *> >::iterator polygon = geometries.begin();
         polygon != geometries.end(); ++polygon) {
        
    }
}

OGRMultiPolygon *PolygonRepair::reconstruct(Triangulation &triangulation) {
    
    // std::cout << "Triangulation: " << triangulation.number_of_faces() << " faces, " << triangulation.number_of_vertices() << " vertices." << std::endl;
    if (triangulation.number_of_faces() < 1) {
        return NULL;
    }
    
    // Reconstruct
    OGRMultiPolygon *outPolygons = new OGRMultiPolygon();
    for (Triangulation::Finite_faces_iterator seedingFace = triangulation.finite_faces_begin(); seedingFace != triangulation.finite_faces_end(); ++seedingFace) {
        
        if (!seedingFace->info().isInInterior()) continue;
        
        // Get boundary
        std::list<Triangulation::Vertex_handle> *vertices = new std::list<Triangulation::Vertex_handle>();
        seedingFace->info().clear();
        if (seedingFace->neighbor(2)->info().isInInterior()) {
            seedingFace->neighbor(2)->info().clear();
            std::list<Triangulation::Vertex_handle> *l2 = getBoundary(seedingFace->neighbor(2), seedingFace->neighbor(2)->index(seedingFace));
            vertices->splice(vertices->end(), *l2);
            delete l2;
        } vertices->push_back(seedingFace->vertex(0));
        if (seedingFace->neighbor(1)->info().isInInterior()) {
            seedingFace->neighbor(1)->info().clear();
            std::list<Triangulation::Vertex_handle> *l1 = getBoundary(seedingFace->neighbor(1), seedingFace->neighbor(1)->index(seedingFace));
            vertices->splice(vertices->end(), *l1);
            delete l1;
        } vertices->push_back(seedingFace->vertex(2));
        if (seedingFace->neighbor(0)->info().isInInterior()) {
            seedingFace->neighbor(0)->info().clear();
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
        outPolygons->addGeometryDirectly(newPolygon);
    } return outPolygons;
}

std::list<Triangulation::Vertex_handle> *PolygonRepair::getBoundary(Triangulation::Face_handle face, int edge) {
    
    std::list<Triangulation::Vertex_handle> *vertices = new std::list<Triangulation::Vertex_handle>();
    
    // Check clockwise edge
    if (!face->is_constrained(face->cw(edge)) && face->neighbor(face->cw(edge))->info().beenReconstructed()) {
		face->neighbor(face->cw(edge))->info().clear();
		std::list<Triangulation::Vertex_handle> *v1 = getBoundary(face->neighbor(face->cw(edge)), face->neighbor(face->cw(edge))->index(face));
		vertices->splice(vertices->end(), *v1);
		delete v1;
	}
	
	// Add central vertex
	vertices->push_back(face->vertex(edge));
	
	// Check counterclockwise edge
    if (!face->is_constrained(face->ccw(edge)) && face->neighbor(face->ccw(edge))->info().beenReconstructed()) {
		face->neighbor(face->ccw(edge))->info().clear();
		std::list<Triangulation::Vertex_handle> *v2 = getBoundary(face->neighbor(face->ccw(edge)), face->neighbor(face->ccw(edge))->index(face));
		vertices->splice(vertices->end(), *v2);
		delete v2;
	}
	
    return vertices;
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