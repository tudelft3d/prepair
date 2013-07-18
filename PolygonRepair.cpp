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
    std::cerr << "PolygonRepair::repairOddEven" << std::endl;
    Triangulation triangulation;
    insertConstraints(triangulation, geometry);
    tagOddEven(triangulation);
    OGRMultiPolygon *outPolygons = reconstruct(triangulation);
    char *outputWKT;
    outPolygons->OGRGeometryCollection::exportToWkt(&outputWKT);
    std::cout << outputWKT << std::endl;
    return outPolygons;
}

OGRMultiPolygon *PolygonRepair::repairPointSet(OGRGeometry *geometry) {
    std::cerr << "PolygonRepair::repairPointSet" << std::endl;
    std::list<std::pair<bool, OGRMultiPolygon *> > repairedRings;   // bool indicates if outer/inner are flipped
    
    switch (geometry->getGeometryType()) {
            
        case wkbLineString: {
            std::cout << "Repairing ring..." << std::endl;
            repairedRings.push_back(std::pair<bool, OGRMultiPolygon *>(false, repairOddEven(geometry)));
            break;
        }
            
        case wkbPolygon: {
            OGRPolygon *polygon = static_cast<OGRPolygon *>(geometry);
            if (polygon->getExteriorRing() != NULL) {
                std::cout << "Repairing polygon outer ring..." << std::endl;
                repairedRings.push_back(std::pair<bool, OGRMultiPolygon *>(false, repairOddEven(polygon->getExteriorRing())));
            } for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing) {
                std::cout << "Repairing polygon inner ring..." << std::endl;
                repairedRings.push_back(std::pair<bool, OGRMultiPolygon *>(true, repairOddEven(polygon->getInteriorRing(currentRing))));
            } break;
        }
            
        case wkbMultiPolygon: {
            OGRMultiPolygon *multipolygon = static_cast<OGRMultiPolygon *>(geometry);
            for (int currentPolygon = 0; currentPolygon < multipolygon->getNumGeometries(); ++currentPolygon) {
                OGRPolygon *polygon = static_cast<OGRPolygon *>(multipolygon->getGeometryRef(currentPolygon));
                if (polygon->getExteriorRing() != NULL) {
                    std::cout << "Repairing multipolygon outer ring..." << std::endl;
                    repairedRings.push_back(std::pair<bool, OGRMultiPolygon *>(false, repairOddEven(polygon->getExteriorRing())));
                } for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing) {
                    std::cout << "Repairing multipolygon inner ring..." << std::endl;
                    repairedRings.push_back(std::pair<bool, OGRMultiPolygon *>(true, repairOddEven(polygon->getInteriorRing(currentRing))));
                }
            } break;
        }
            
        default:
            std::cerr << "PolygonRepair::repairPointSet: Cannot understand input. Only polygons are supported." << std::endl;
            return new OGRMultiPolygon();
            break;
    }
    
    Triangulation triangulation;
    for (std::list<std::pair<bool, OGRMultiPolygon *> >::iterator currentMultipolygon = repairedRings.begin();
         currentMultipolygon != repairedRings.end(); ++currentMultipolygon)
        insertConstraints(triangulation, currentMultipolygon->second);
    printEdges(triangulation);
    tagPointSet(triangulation, repairedRings);
    OGRMultiPolygon *outPolygons = reconstruct(triangulation);
    return outPolygons;
}

void PolygonRepair::removeSmallPolygons(OGRMultiPolygon *outPolygons, double minArea) {
    std::cout << "Removing polygons smaller than " << minArea << " unit^2." << std::endl;
    for (int currentPolygon = outPolygons->getNumGeometries()-1; currentPolygon >= 0; --currentPolygon) {
        if (static_cast<OGRPolygon *>(outPolygons->getGeometryRef(currentPolygon))->get_Area() < minArea) {
            outPolygons->removeGeometry(currentPolygon);
        }
    }
}

OGRMultiPolygon *PolygonRepair::isr(OGRGeometry *geometry) {
    
}

void PolygonRepair::insertConstraints(Triangulation &triangulation, OGRGeometry* geometry) {
    Triangulation::Vertex_handle va, vb;
    Triangulation::Face_handle faceOfEdge;
    int indexOfEdge;
    
    switch (geometry->getGeometryType()) {
            
        case wkbLineString: {
            OGRLinearRing *ring = static_cast<OGRLinearRing *>(geometry);
            
            for (int currentPoint = 0; currentPoint < ring->getNumPoints(); ++currentPoint) {
                va = triangulation.insert(Point(ring->getX(currentPoint),
                                                ring->getY(currentPoint)));
                vb = triangulation.insert(Point(ring->getX((currentPoint+1)%ring->getNumPoints()),
                                                ring->getY((currentPoint+1)%ring->getNumPoints())));
                std::cerr << "Inserting contraint <" << va->point() << ", " << vb->point() << "> in outer ring..." << std::endl;
                if (triangulation.is_edge(va, vb, faceOfEdge, indexOfEdge)) {
                    if (triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(faceOfEdge, indexOfEdge))) {
                        std::cerr << "\tRemoving duplicate constraint <" << va->point() << ", " << vb->point() << "> in outer ring..." << std::endl;
                        triangulation.insert_constraint(va, vb);
                        triangulation.remove_constraint(va, vb);
                    } else triangulation.insert_constraint(va, vb);
                } else triangulation.insert_constraint(va, vb);
            } break;
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
                        std::cerr << "\tRemoving duplicate constraint <" << va->point() << ", " << vb->point() << "> in outer ring..." << std::endl;
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
                            std::cerr << "\tRemoving duplicate constraint <" << va->point() << ", " << vb->point() << "> in inner ring..." << std::endl;
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
                            std::cerr << "\tRemoving duplicate constraint <" << va->point() << ", " << vb->point() << "> in outer ring..." << std::endl;
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
                                std::cerr << "\tRemoving duplicate constraint <" << va->point() << ", " << vb->point() << "> in inner ring..." << std::endl;
                                triangulation.insert_constraint(va, vb);
                                triangulation.remove_constraint(va, vb);
                            } else triangulation.insert_constraint(va, vb);
                        } else triangulation.insert_constraint(va, vb);
                    }
                }
            } break;
        }
        
        default:
            std::cerr << "PolygonRepair::insertConstraints: Cannot understand input. Only polygons are supported." << std::endl;
            return;
            break;
    }
    
    // Remove partially even-overlapping subconstraints
    std::cerr << "Checking subconstraints..." << std::endl;
    for (Triangulation::Subconstraint_iterator currentEdge = triangulation.subconstraints_begin();
         currentEdge != triangulation.subconstraints_end();
         ++currentEdge) {
        if (triangulation.number_of_enclosing_constraints(currentEdge->first.first, currentEdge->first.second) % 2 == 0) {
            std::cerr << "\tRemoving even-overlapped subconstraint <" << currentEdge->first.first->point() << ", " << currentEdge->first.second->point() << ">..." << std::endl;
            triangulation.remove_constraint(currentEdge->first.first, currentEdge->first.second);
        }
    }
}

void PolygonRepair::tagOddEven(Triangulation &triangulation) {
    std::cout << "PolygonRepair::tagOddEven" << std::endl;
	
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
            if (labellingInterior) {
                std::cout << "Tagged ";
                printTriangle(currentFace);
                std::cout << std::endl;
            } for (int currentEdge = 0; currentEdge < 3; ++currentEdge) {
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
    std::cout << "PolygonRepair::tagPointSet" << std::endl;
    
    // Clean tags
    for (Triangulation::Face_handle currentFace = triangulation.all_faces_begin(); currentFace != triangulation.all_faces_end(); ++currentFace)
        currentFace->info().clear();
    
    std::stack<Triangulation::Face_handle> borderTriangles, taggingStack, untaggingStack;
    Triangulation::Vertices_in_constraint_iterator currentVertex, nextVertex, lastVertex;
    Triangulation::Vertex_handle va, vb;
    Triangulation::Face_handle faceOfEdge, faceOfSubedge;
    int indexOfEdge, indexOfSubedge;
    bool sameOrder;
    
    // Add all repaired outer rings
    std::cout << "Flooding outer rings..." << std::endl;
    for (std::list<std::pair<bool, OGRMultiPolygon *> >::iterator multipolygon = geometries.begin();
         multipolygon != geometries.end(); ++multipolygon) {
        for (int currentPolygon = 0; currentPolygon < multipolygon->second->getNumGeometries(); ++currentPolygon) {
            OGRPolygon *polygon = static_cast<OGRPolygon *>(multipolygon->second->getGeometryRef(currentPolygon));
            taggingStack = std::stack<Triangulation::Face_handle>();
            
            // Outer
            if (!multipolygon->first) for (int currentPoint = 0; currentPoint < polygon->getExteriorRing()->getNumPoints(); ++currentPoint) {
                va = triangulation.insert(Point(polygon->getExteriorRing()->getX(currentPoint),
                                                polygon->getExteriorRing()->getY(currentPoint)));
                vb = triangulation.insert(Point(polygon->getExteriorRing()->getX((currentPoint+1)%polygon->getExteriorRing()->getNumPoints()),
                                                polygon->getExteriorRing()->getY((currentPoint+1)%polygon->getExteriorRing()->getNumPoints())));
                std::cerr << "Tagging contraint <" << va->point() << ", " << vb->point() << "> in outer ring..." << std::endl;
                if (triangulation.is_edge(va, vb, faceOfEdge, indexOfEdge)) {
                    if (triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(faceOfEdge, indexOfEdge))) {
                        currentVertex = triangulation.vertices_in_constraint_begin(va, vb);
                        nextVertex = currentVertex;
                        ++nextVertex;
                        lastVertex = triangulation.vertices_in_constraint_end(va, vb);
                        if (*currentVertex == va) sameOrder = true;
                        else sameOrder = false;
                        while (nextVertex != lastVertex) {
                            if (!sameOrder) {
                                if (!triangulation.is_edge(*currentVertex, *nextVertex, faceOfSubedge, indexOfSubedge)) {
                                    std::cerr << "PolygonRepair::tagPointSet: Cannot find edge!" << std::endl;
                                    return;
                                }
                            } else {
                                if (!triangulation.is_edge(*nextVertex, *currentVertex, faceOfSubedge, indexOfSubedge)) {
                                    std::cerr << "PolygonRepair::tagPointSet: Cannot find edge!" << std::endl;
                                    return;
                                }
                            } currentVertex = nextVertex;
                            ++nextVertex;
                            std::cout << "Marked in internal border: ";
                            printTriangle(faceOfSubedge);
                            std::cout << std::endl;
                            std::cout << "Marked in external border: ";
                            printTriangle(faceOfSubedge->neighbor(indexOfSubedge));
                            std::cout << std::endl;
                            borderTriangles.push(faceOfSubedge);
                            faceOfSubedge->info().isOnBorder(true);
                            borderTriangles.push(faceOfSubedge->neighbor(indexOfSubedge));
                            faceOfSubedge->neighbor(indexOfSubedge)->info().isOnBorder(true);
                            taggingStack.push(faceOfSubedge);
                        }
                    } else {
                        std::cout << "\tNot constrained anymore. Skipping..." << std::endl;
                    }
                } else {
                    std::cout << "\tNot an edge anymore. Skipping..." << std::endl;
                }
            }
            
            // Inner
            else for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing) {
                for (int currentPoint = 0; currentPoint < polygon->getInteriorRing(currentRing)->getNumPoints(); ++currentPoint) {
                    va = triangulation.insert(Point(polygon->getInteriorRing(currentRing)->getX(currentPoint),
                                                    polygon->getInteriorRing(currentRing)->getY(currentPoint)));
                    vb = triangulation.insert(Point(polygon->getInteriorRing(currentRing)->getX((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints()),
                                                    polygon->getInteriorRing(currentRing)->getY((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints())));
                    std::cerr << "Tagging contraint <" << va->point() << ", " << vb->point() << "> in inner ring..." << std::endl;
                    if (triangulation.is_edge(va, vb, faceOfEdge, indexOfEdge)) {
                        if (triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(faceOfEdge, indexOfEdge))) {
                            currentVertex = triangulation.vertices_in_constraint_begin(va, vb);
                            nextVertex = currentVertex;
                            ++nextVertex;
                            lastVertex = triangulation.vertices_in_constraint_end(va, vb);
                            if (*currentVertex == va) sameOrder = true;
                            else sameOrder = false;
                            while (nextVertex != lastVertex) {
                                if (sameOrder) {
                                    if (!triangulation.is_edge(*currentVertex, *nextVertex, faceOfSubedge, indexOfSubedge)) {
                                        std::cerr << "PolygonRepair::tagPointSet: Cannot find edge!" << std::endl;
                                        return;
                                    }
                                } else {
                                    if (!triangulation.is_edge(*nextVertex, *currentVertex, faceOfSubedge, indexOfSubedge)) {
                                        std::cerr << "PolygonRepair::tagPointSet: Cannot find edge!" << std::endl;
                                        return;
                                    }
                                } currentVertex = nextVertex;
                                ++nextVertex;
                                std::cout << "Marked in internal border: ";
                                printTriangle(faceOfSubedge);
                                std::cout << std::endl;
                                std::cout << "Marked in external border: ";
                                printTriangle(faceOfSubedge->neighbor(indexOfSubedge));
                                std::cout << std::endl;
                                borderTriangles.push(faceOfSubedge);
                                faceOfSubedge->info().isOnBorder(true);
                                borderTriangles.push(faceOfSubedge->neighbor(indexOfSubedge));
                                faceOfSubedge->neighbor(indexOfSubedge)->info().isOnBorder(true);
                                taggingStack.push(faceOfSubedge);
                            }
                        } else {
                            std::cout << "\tNot constrained anymore. Skipping..." << std::endl;
                        }
                    } else {
                        std::cout << "\tNot an edge anymore. Skipping..." << std::endl;
                    }
                }
            }
            
            // Expand the tags
            untaggingStack = std::stack<Triangulation::Face_handle>();
            while (!taggingStack.empty()) {
                Triangulation::Face_handle currentFace = taggingStack.top();
                taggingStack.pop();
                if (currentFace->info().beenTagged()) continue;
                std::cout << "Tagging ";
                printTriangle(currentFace);
                std::cout << std::endl;
                currentFace->info().isInInterior(true);
                untaggingStack.push(currentFace);
                for (int currentEdge = 0; currentEdge < 3; ++currentEdge) {
                    if (!currentFace->neighbor(currentEdge)->info().isOnBorder() &&
                        !currentFace->neighbor(currentEdge)->info().beenTagged()) {
                        taggingStack.push(currentFace->neighbor(currentEdge));
                    }
                }
            }
            
            // Remove border tags
            while (!borderTriangles.empty()) {
                borderTriangles.top()->info().isOnBorder(false);
                borderTriangles.pop();
            }
            
            // Remove tagged tags
            while (!untaggingStack.empty()) {
                untaggingStack.top()->info().beenTagged(false);
                untaggingStack.pop();
            }
        }
    }
    
    // TODO: Subtract all repaired inner rings
    std::cout << "Carving inner rings..." << std::endl;
    for (std::list<std::pair<bool, OGRMultiPolygon *> >::iterator multipolygon = geometries.begin();
         multipolygon != geometries.end(); ++multipolygon) {
        for (int currentPolygon = 0; currentPolygon < multipolygon->second->getNumGeometries(); ++currentPolygon) {
            OGRPolygon *polygon = static_cast<OGRPolygon *>(multipolygon->second->getGeometryRef(currentPolygon));
            taggingStack = std::stack<Triangulation::Face_handle>();
            
            // Outer
            if (multipolygon->first) for (int currentPoint = 0; currentPoint < polygon->getExteriorRing()->getNumPoints(); ++currentPoint) {
                va = triangulation.insert(Point(polygon->getExteriorRing()->getX(currentPoint),
                                                polygon->getExteriorRing()->getY(currentPoint)));
                vb = triangulation.insert(Point(polygon->getExteriorRing()->getX((currentPoint+1)%polygon->getExteriorRing()->getNumPoints()),
                                                polygon->getExteriorRing()->getY((currentPoint+1)%polygon->getExteriorRing()->getNumPoints())));
                std::cerr << "Carving contraint <" << va->point() << ", " << vb->point() << "> in outer ring..." << std::endl;
                if (triangulation.is_edge(va, vb, faceOfEdge, indexOfEdge)) {
                    if (triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(faceOfEdge, indexOfEdge))) {
                        currentVertex = triangulation.vertices_in_constraint_begin(va, vb);
                        nextVertex = currentVertex;
                        ++nextVertex;
                        lastVertex = triangulation.vertices_in_constraint_end(va, vb);
                        if (*currentVertex == va) sameOrder = true;
                        else sameOrder = false;
                        while (nextVertex != lastVertex) {
                            if (!sameOrder) {
                                if (!triangulation.is_edge(*currentVertex, *nextVertex, faceOfSubedge, indexOfSubedge)) {
                                    std::cerr << "PolygonRepair::tagPointSet: Cannot find edge!" << std::endl;
                                    return;
                                }
                            } else {
                                if (!triangulation.is_edge(*nextVertex, *currentVertex, faceOfSubedge, indexOfSubedge)) {
                                    std::cerr << "PolygonRepair::tagPointSet: Cannot find edge!" << std::endl;
                                    return;
                                }
                            } currentVertex = nextVertex;
                            ++nextVertex;
                            std::cout << "Marked in internal border: ";
                            printTriangle(faceOfSubedge);
                            std::cout << std::endl;
                            std::cout << "Marked in external border: ";
                            printTriangle(faceOfSubedge->neighbor(indexOfSubedge));
                            std::cout << std::endl;
                            borderTriangles.push(faceOfSubedge);
                            faceOfSubedge->info().isOnBorder(true);
                            borderTriangles.push(faceOfSubedge->neighbor(indexOfSubedge));
                            faceOfSubedge->neighbor(indexOfSubedge)->info().isOnBorder(true);
                            taggingStack.push(faceOfSubedge);
                        }
                    } else {
                        std::cout << "\tNot constrained anymore. Skipping..." << std::endl;
                    }
                } else {
                    std::cout << "\tNot an edge anymore. Skipping..." << std::endl;
                }
            }
            
            // Inner
            else for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing) {
                for (int currentPoint = 0; currentPoint < polygon->getInteriorRing(currentRing)->getNumPoints(); ++currentPoint) {
                    va = triangulation.insert(Point(polygon->getInteriorRing(currentRing)->getX(currentPoint),
                                                    polygon->getInteriorRing(currentRing)->getY(currentPoint)));
                    vb = triangulation.insert(Point(polygon->getInteriorRing(currentRing)->getX((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints()),
                                                    polygon->getInteriorRing(currentRing)->getY((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints())));
                    std::cerr << "Tagging contraint <" << va->point() << ", " << vb->point() << "> in inner ring..." << std::endl;
                    if (triangulation.is_edge(va, vb, faceOfEdge, indexOfEdge)) {
                        if (triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(faceOfEdge, indexOfEdge))) {
                            currentVertex = triangulation.vertices_in_constraint_begin(va, vb);
                            nextVertex = currentVertex;
                            ++nextVertex;
                            lastVertex = triangulation.vertices_in_constraint_end(va, vb);
                            if (*currentVertex == va) sameOrder = true;
                            else sameOrder = false;
                            while (nextVertex != lastVertex) {
                                if (sameOrder) {
                                    if (!triangulation.is_edge(*currentVertex, *nextVertex, faceOfSubedge, indexOfSubedge)) {
                                        std::cerr << "PolygonRepair::tagPointSet: Cannot find edge!" << std::endl;
                                        return;
                                    }
                                } else {
                                    if (!triangulation.is_edge(*nextVertex, *currentVertex, faceOfSubedge, indexOfSubedge)) {
                                        std::cerr << "PolygonRepair::tagPointSet: Cannot find edge!" << std::endl;
                                        return;
                                    }
                                } currentVertex = nextVertex;
                                ++nextVertex;
                                std::cout << "Marked in internal border: ";
                                printTriangle(faceOfSubedge);
                                std::cout << std::endl;
                                std::cout << "Marked in external border: ";
                                printTriangle(faceOfSubedge->neighbor(indexOfSubedge));
                                std::cout << std::endl;
                                borderTriangles.push(faceOfSubedge);
                                faceOfSubedge->info().isOnBorder(true);
                                borderTriangles.push(faceOfSubedge->neighbor(indexOfSubedge));
                                faceOfSubedge->neighbor(indexOfSubedge)->info().isOnBorder(true);
                                taggingStack.push(faceOfSubedge);
                            }
                        } else {
                            std::cout << "\tNot constrained anymore. Skipping..." << std::endl;
                        }
                    } else {
                        std::cout << "\tNot an edge anymore. Skipping..." << std::endl;
                    }
                }
            }
            
            // Expand the tags
            untaggingStack = std::stack<Triangulation::Face_handle>();
            while (!taggingStack.empty()) {
                Triangulation::Face_handle currentFace = taggingStack.top();
                taggingStack.pop();
                if (currentFace->info().beenTagged()) continue;
                currentFace->info().isInInterior(false);
                untaggingStack.push(currentFace);
                for (int currentEdge = 0; currentEdge < 3; ++currentEdge) {
                    if (!currentFace->neighbor(currentEdge)->info().isOnBorder() &&
                        !currentFace->neighbor(currentEdge)->info().beenTagged()) {
                        taggingStack.push(currentFace->neighbor(currentEdge));
                    }
                }
            }
            
            // Remove border tags
            while (!borderTriangles.empty()) {
                borderTriangles.top()->info().isOnBorder(false);
                borderTriangles.pop();
            }
            
            // Remove tagged tags
            while (!untaggingStack.empty()) {
                untaggingStack.top()->info().beenTagged(false);
                untaggingStack.pop();
            }
        }
    }
}

OGRMultiPolygon *PolygonRepair::reconstruct(Triangulation &triangulation) {
    std::cout << "PolygonRepair::reconstruct" << std::endl;
    
    // std::cout << "Triangulation: " << triangulation.number_of_faces() << " faces, " << triangulation.number_of_vertices() << " vertices." << std::endl;
    if (triangulation.number_of_faces() < 1) {
        return NULL;
    }
    
    // Reconstruct
    OGRMultiPolygon *outPolygons = new OGRMultiPolygon();
    for (Triangulation::Finite_faces_iterator seedingFace = triangulation.finite_faces_begin(); seedingFace != triangulation.finite_faces_end(); ++seedingFace) {
        
        if (!seedingFace->info().isInInterior() || seedingFace->info().beenReconstructed()) continue;
        std::cout << "Seeding face: ";
        printTriangle(seedingFace);
        std::cout << std::endl;
        seedingFace->info().beenReconstructed(true);
        if (!seedingFace->info().beenReconstructed()) {
            std::cout << "ERROR! Should be marked as reconstructed!!!" << std::endl;
        }
        
        // Get boundary
        std::list<Triangulation::Vertex_handle> *vertices = new std::list<Triangulation::Vertex_handle>();
        if (seedingFace->neighbor(2)->info().isInInterior() && !seedingFace->neighbor(2)->info().beenReconstructed()) {
            seedingFace->neighbor(2)->info().beenReconstructed(true);
            std::list<Triangulation::Vertex_handle> *l2 = getBoundary(seedingFace->neighbor(2), seedingFace->neighbor(2)->index(seedingFace));
            std::cout << "Putting ";
            printChain(*l2);
            vertices->splice(vertices->end(), *l2);
            delete l2;
        } std::cout << "Inserting " << seedingFace->vertex(0)->point() << std::endl;
        vertices->push_back(seedingFace->vertex(0));
        if (seedingFace->neighbor(1)->info().isInInterior() && !seedingFace->neighbor(1)->info().beenReconstructed()) {
            seedingFace->neighbor(1)->info().beenReconstructed(true);
            std::list<Triangulation::Vertex_handle> *l1 = getBoundary(seedingFace->neighbor(1), seedingFace->neighbor(1)->index(seedingFace));
            std::cout << "Putting ";
            printChain(*l1);
            vertices->splice(vertices->end(), *l1);
            delete l1;
        } std::cout << "Inserting " << seedingFace->vertex(2)->point() << std::endl;
        vertices->push_back(seedingFace->vertex(2));
        if (seedingFace->neighbor(0)->info().isInInterior() && !seedingFace->neighbor(0)->info().beenReconstructed()) {
            seedingFace->neighbor(0)->info().beenReconstructed(true);
            std::list<Triangulation::Vertex_handle> *l0 = getBoundary(seedingFace->neighbor(0), seedingFace->neighbor(0)->index(seedingFace));
            std::cout << "Putting ";
            printChain(*l0);
            vertices->splice(vertices->end(), *l0);
            delete l0;
        } std::cout << "Inserting " << seedingFace->vertex(1)->point() << std::endl;
        vertices->push_back(seedingFace->vertex(1));
        std::cout << "Chain: ";
        printChain(*vertices);
        
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
    std::cout << "PolygonRepair::getBoundary" << std::endl;
    
    std::list<Triangulation::Vertex_handle> *vertices = new std::list<Triangulation::Vertex_handle>();
    
    // Check clockwise edge
    //std::cout << "CW" << std::endl;
    if (face->neighbor(face->cw(edge))->info().isInInterior() && !face->neighbor(face->cw(edge))->info().beenReconstructed()) {
		face->neighbor(face->cw(edge))->info().beenReconstructed(true);
		std::list<Triangulation::Vertex_handle> *v1 = getBoundary(face->neighbor(face->cw(edge)), face->neighbor(face->cw(edge))->index(face));
		vertices->splice(vertices->end(), *v1);
		delete v1;
	}
	
	// Add central vertex
    //std::cout << "CV: " << face->vertex(edge)->point() << std::endl;
	vertices->push_back(face->vertex(edge));
	
	// Check counterclockwise edge
    //std::cout << "CCW" << std::endl;
    if (face->neighbor(face->ccw(edge))->info().isInInterior() && !face->neighbor(face->ccw(edge))->info().beenReconstructed()) {
		face->neighbor(face->ccw(edge))->info().beenReconstructed(true);
		std::list<Triangulation::Vertex_handle> *v2 = getBoundary(face->neighbor(face->ccw(edge)), face->neighbor(face->ccw(edge))->index(face));
		vertices->splice(vertices->end(), *v2);
		delete v2;
	}
	
    //std::cout << "Done." << std::endl;
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

void PolygonRepair::printTriangle(Triangulation::Face_handle triangle) {
    std::cout << "TRIANGLE(" << triangle->vertex(0)->point() << ", " << triangle->vertex(1)->point() << ", " << triangle->vertex(2)->point() << ")";
}

void PolygonRepair::printChain(std::list<Triangulation::Vertex_handle> &chain) {
    for (std::list<Triangulation::Vertex_handle>::iterator currentVertex = chain.begin(); currentVertex != chain.end(); ++currentVertex) {
        std::cout << (*currentVertex)->point() << ", ";
    } std::cout << std::endl;
}