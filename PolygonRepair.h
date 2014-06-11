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

#include "definitions.h"

class PolygonRepair {
public:
  void repairOddEven(MultiPolygon &inGeometry, MultiPolygon &outGeometry, bool timeResults = false);
//  void repairPointSet(MultiPolygon &inGeometry, MultiPolygon &outGeometry, bool timeResults = false);
//  void removeSmallPolygons(MultiPolygon &geometry, K::FT minArea);
//  void isr(MultiPolygon &inGeometry, MultiPolygon &outGeometry, K::FT tolerance);
//  K::FT computeRobustness(MultiPolygon &geometry);
//  bool saveToShp(OGRGeometry* geometry, const char *fileName);
  
private:
  Triangulation triangulation;
  
  void insertConstraints(Triangulation &triangulation, MultiPolygon &geometry, bool removeOverlappingConstraints = true);
  void tagOddEven(Triangulation &triangulation);
//  void tagPointSet(Triangulation &triangulation, std::list<std::pair<bool, OGRMultiPolygon *> > &geometries);
  void reconstruct(Triangulation &triangulation, MultiPolygon &outPolygons);
  void getBoundary(Triangulation::Face_handle face, int edge, LinearRing &vertices);
  
  // Debug functions
  void printEdges(Triangulation &triangulation);
  void printTriangle(Triangulation::Face_handle triangle);
  void printChain(LinearRing &chain);
};

#endif