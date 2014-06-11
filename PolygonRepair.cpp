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

#include "PolygonRepair.h"

void PolygonRepair::ogrToMultiPolygon(OGRGeometry *inGeometry, MultiPolygon &outGeometry) {
  outGeometry.clear();
  switch (inGeometry->getGeometryType()) {
      
    case wkbLineString: {
      OGRLinearRing *geometry = static_cast<OGRLinearRing *>(inGeometry);
      outGeometry.push_back(Polygon());
      MultiPolygon::reverse_iterator newPolygon = outGeometry.rbegin();
      newPolygon->push_back(LinearRing());
      Polygon::reverse_iterator newRing = newPolygon->rbegin();
      for (int currentVertex = 0; currentVertex < geometry->getNumPoints(); ++currentVertex) {
        newRing->push
      }
    }
  }
}

void PolygonRepair::repairOddEven(MultiPolygon &inGeometry, MultiPolygon &outGeometry, bool timeResults) {
  triangulation.clear();
  time_t thisTime, totalTime;
  thisTime = time(NULL);
  insertConstraints(triangulation, inGeometry);
  totalTime = time(NULL)-thisTime;
  if (timeResults) std::cout << "Triangulation: " << totalTime/60 << " minutes " << totalTime%60 << " seconds." << std::endl;
  thisTime = time(NULL);
  tagOddEven(triangulation);
  totalTime = time(NULL)-thisTime;
  if (timeResults) std::cout << "Tagging: " << totalTime/60 << " minutes " << totalTime%60 << " seconds." << std::endl;
  thisTime = time(NULL);
  reconstruct(triangulation, outGeometry);
  totalTime = time(NULL)-thisTime;
  if (timeResults) std::cout << "Reconstruction: " << totalTime/60 << " minutes " << totalTime%60 << " seconds." << std::endl;
}

//void PolygonRepair::repairPointSet(OGRGeometry *geometry, MultiPolygon &outPolygons, bool timeResults) {
//  triangulation.clear();
//  std::list<MultiPolygon> repairedRings, flippedRepairedRings;
//  
//  time_t thisTime, totalTime;
//  thisTime = time(NULL);
//  switch (geometry->getGeometryType()) {
//      
//    case wkbLineString: {
//      repairedRings.push_back(MultiPolygon());
//      repairOddEven(geometry, repairedRings.back(), false);
//      break;
//    }
//
//    case wkbPolygon: {
//      OGRPolygon *polygon = static_cast<OGRPolygon *>(geometry);
//      if (polygon->getExteriorRing() != NULL) {
//        repairedRings.push_back(MultiPolygon());
//        repairOddEven(polygon->getExteriorRing(), repairedRings.back(), false);
//      } for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing) {
//        flippedRepairedRings.push_back(MultiPolygon());
//        repairOddEven(polygon->getInteriorRing(currentRing), flippedRepairedRings.back(), false);
//      } break;
//    }
//      
//    case wkbMultiPolygon: {
//      OGRMultiPolygon *multipolygon = static_cast<OGRMultiPolygon *>(geometry);
//      for (int currentPolygon = 0; currentPolygon < multipolygon->getNumGeometries(); ++currentPolygon) {
//        OGRPolygon *polygon = static_cast<OGRPolygon *>(multipolygon->getGeometryRef(currentPolygon));
//        if (polygon->getExteriorRing() != NULL) {
//          repairedRings.push_back(MultiPolygon());
//          repairOddEven(polygon->getExteriorRing(), repairedRings.back(), false);
//        } for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing) {
//          flippedRepairedRings.push_back(MultiPolygon());
//          repairOddEven(polygon->getInteriorRing(currentRing), flippedRepairedRings.back(), false);
//        }
//      } break;
//    }
//      
//    default:
//      std::cerr << "PolygonRepair::repairPointSet: Cannot understand input." << std::endl;
//      return;
//  } totalTime = time(NULL)-thisTime;
//  if (timeResults) std::cout << "Repairing individual rings: " << totalTime/60 << " minutes " << totalTime%60 << " seconds." << std::endl;
//  
//  thisTime = time(NULL);
//  for (std::list<MultiPolygon>::iterator currentMultipolygon = repairedRings.begin();
//       currentMultipolygon != repairedRings.end(); ++currentMultipolygon)
//    insertConstraints(triangulation, currentMultipolygon, false);
//  for (std::list<MultiPolygon>::iterator currentMultipolygon = flippedRepairedRings.begin();
//       currentMultipolygon != flippedRepairedRings.end(); ++currentMultipolygon)
//    insertConstraints(triangulation, currentMultipolygon, false);
//  totalTime = time(NULL)-thisTime;
//  if (timeResults) std::cout << "Triangulation: " << totalTime/60 << " minutes " << totalTime%60 << " seconds." << std::endl;
//  thisTime = time(NULL);
//  tagPointSet(triangulation, repairedRings, flippedRepairedRings);
//  totalTime = time(NULL)-thisTime;
//  if (timeResults) std::cout << "Tagging: " << totalTime/60 << " minutes " << totalTime%60 << " seconds." << std::endl;
//  thisTime = time(NULL);
//  reconstruct(triangulation, outPolygons);
//  totalTime = time(NULL)-thisTime;
//  if (timeResults) std::cout << "Reconstruction: " << totalTime/60 << " minutes " << totalTime%60 << " seconds." << std::endl;
//  return outPolygons;
//}

//void PolygonRepair::removeSmallPolygons(MultiPolygon &outPolygons, double minArea) {
//  //std::cout << "Removing polygons smaller than " << minArea << " unit^2." << std::endl;
//  for (int currentPolygon = outPolygons->getNumGeometries()-1; currentPolygon >= 0; --currentPolygon) {
//    if (static_cast<OGRPolygon *>(outPolygons->getGeometryRef(currentPolygon))->get_Area() < minArea) {
//      outPolygons->removeGeometry(currentPolygon);
//    }
//  }
//}

//OGRMultiPolygon *PolygonRepair::isr(OGRGeometry *geometry, double tolerance) {
//  
//  // Convert to multipolygon
//  OGRMultiPolygon *multipolygon = new OGRMultiPolygon();
//  switch (geometry->getGeometryType()) {
//    case wkbLineString: {
//      OGRLineString *oldRing = static_cast<OGRLineString *>(geometry);
//      OGRPolygon *polygon = new OGRPolygon();
//      OGRLinearRing *ring = new OGRLinearRing();
//      OGRPoint lastPoint, point;
//      for (int currentPoint = 0; currentPoint < oldRing->getNumPoints(); ++currentPoint) {
//        oldRing->getPoint(currentPoint, &point);
//        if (currentPoint > 0) {
//          if (!point.Equals(&lastPoint)) ring->addPoint(&point);
//        } lastPoint = point;
//      } polygon->addRingDirectly(ring);
//      multipolygon->addGeometryDirectly(polygon);
//      break;
//    }
//      
//    case wkbPolygon: {
//      OGRPolygon *oldPolygon = static_cast<OGRPolygon *>(geometry);
//      OGRPolygon *polygon = new OGRPolygon();
//      OGRLinearRing *oldRing = oldPolygon->getExteriorRing();
//      OGRLinearRing *ring = new OGRLinearRing();
//      OGRPoint lastPoint, point;
//      for (int currentPoint = 0; currentPoint < oldRing->getNumPoints(); ++currentPoint) {
//        oldRing->getPoint(currentPoint, &point);
//        if (currentPoint > 0) {
//          if (!point.Equals(&lastPoint)) ring->addPoint(&point);
//        } lastPoint = point;
//      } polygon->addRingDirectly(ring);
//      for (int currentRing = 0; currentRing < oldPolygon->getNumInteriorRings(); ++currentRing) {
//        oldRing = oldPolygon->getInteriorRing(currentRing);
//        ring = new OGRLinearRing();
//        for (int currentPoint = 0; currentPoint < oldRing->getNumPoints(); ++currentPoint) {
//          oldRing->getPoint(currentPoint, &point);
//          if (currentPoint > 0) {
//            if (!point.Equals(&lastPoint)) ring->addPoint(&point);
//          } lastPoint = point;
//        } polygon->addRingDirectly(ring);
//      } multipolygon->addGeometryDirectly(polygon);
//      break;
//    }
//      
//    case wkbMultiPolygon: {
//      OGRMultiPolygon *oldMultipolygon = static_cast<OGRMultiPolygon *>(geometry);
//      OGRPolygon *oldPolygon, *polygon;
//      OGRLinearRing *oldRing, *ring;
//      OGRPoint lastPoint, point;
//      for (int currentPolygon = 0; currentPolygon < oldMultipolygon->getNumGeometries(); ++currentPolygon) {
//        oldPolygon = static_cast<OGRPolygon *>(oldMultipolygon->getGeometryRef(currentPolygon));
//        polygon = new OGRPolygon();
//        oldRing = oldPolygon->getExteriorRing();
//        ring = new OGRLinearRing();
//        for (int currentPoint = 0; currentPoint < oldRing->getNumPoints(); ++currentPoint) {
//          oldRing->getPoint(currentPoint, &point);
//          if (currentPoint > 0) {
//            if (!point.Equals(&lastPoint)) ring->addPoint(&point);
//          } lastPoint = point;
//        } polygon->addRingDirectly(ring);
//        for (int currentRing = 0; currentRing < oldPolygon->getNumInteriorRings(); ++currentRing) {
//          oldRing = oldPolygon->getInteriorRing(currentRing);
//          ring = new OGRLinearRing();
//          for (int currentPoint = 0; currentPoint < oldRing->getNumPoints(); ++currentPoint) {
//            oldRing->getPoint(currentPoint, &point);
//            if (currentPoint > 0) {
//              if (!point.Equals(&lastPoint)) ring->addPoint(&point);
//            } lastPoint = point;
//          } polygon->addRingDirectly(ring);
//        } multipolygon->addGeometryDirectly(polygon);
//      } break;
//    }
//      
//    default:
//      std::cerr << "PolygonRepair::isr: Cannot understand input." << std::endl;
//      return new OGRMultiPolygon();
//      break;
//  } multipolygon->closeRings();
//  
//  // Put all line segments in a list (in order)
//  // TODO: Only do point location once for segment.
//  std::list<Segment> inputSegments;
//  for (int currentPolygon = 0; currentPolygon < multipolygon->getNumGeometries(); ++currentPolygon) {
//    OGRPolygon *polygon = static_cast<OGRPolygon *>(multipolygon->getGeometryRef(currentPolygon));
//    for (int currentPoint = 1; currentPoint < polygon->getExteriorRing()->getNumPoints(); ++currentPoint) {
//      Segment thisSegment(Point(polygon->getExteriorRing()->getX(currentPoint-1), polygon->getExteriorRing()->getY(currentPoint-1)),
//                          Point(polygon->getExteriorRing()->getX(currentPoint), polygon->getExteriorRing()->getY(currentPoint)));
//      if (thisSegment.is_degenerate()) {
//        std::cout << "This should not happen!!!" << std::endl;
//        return new OGRMultiPolygon();
//      } inputSegments.push_back(thisSegment);
//    } for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing) {
//      for (int currentPoint = 1; currentPoint < polygon->getInteriorRing(currentRing)->getNumPoints(); ++currentPoint) {
//        Segment thisSegment(Point(polygon->getExteriorRing()->getX(currentPoint-1), polygon->getExteriorRing()->getY(currentPoint-1)),
//                            Point(polygon->getExteriorRing()->getX(currentPoint), polygon->getExteriorRing()->getY(currentPoint)));
//        if (thisSegment.is_degenerate()) {
//          std::cout << "This should not happen!!!" << std::endl;
//          return new OGRMultiPolygon();
//        } std::cout << thisSegment << std::endl;
//        inputSegments.push_back(thisSegment);
//      }
//    }
//  }
//  
//  // Snap
//  std::list<std::list<Point> > outputPolylines;
//  CGAL::snap_rounding_2<SRT, std::list<Segment>::const_iterator, std::list<std::list<Point> > >(inputSegments.begin(), inputSegments.end(), outputPolylines, tolerance, true, false, 2);
//  
//  // Put it back into an OGRMultiPolygon
//  OGRMultiPolygon *snappedMultipolygon = new OGRMultiPolygon();
//  std::list<std::list<Point> >::iterator currentPolyline = outputPolylines.begin();
//  for (int currentPolygon = 0; currentPolygon < multipolygon->getNumGeometries(); ++currentPolygon) {
//    OGRPolygon *polygon = static_cast<OGRPolygon *>(multipolygon->getGeometryRef(currentPolygon));
//    OGRPolygon *snappedPolygon = new OGRPolygon();
//    OGRLinearRing *snappedOuterRing = new OGRLinearRing();
//    for (int currentPoint = 1; currentPoint < polygon->getExteriorRing()->getNumPoints(); ++currentPoint) {
//      std::list<Point>::iterator point = currentPolyline->begin();
//      ++point;
//      while (point != currentPolyline->end()) {
//        snappedOuterRing->addPoint(CGAL::to_double(point->x()), CGAL::to_double(point->y()));
//        ++point;
//      } ++currentPolyline;
//    } snappedPolygon->addRingDirectly(snappedOuterRing);
//    for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing) {
//      OGRLinearRing *snappedInnerRing = new OGRLinearRing();
//      for (int currentPoint = 1; currentPoint < polygon->getInteriorRing(currentRing)->getNumPoints(); ++currentPoint) {
//        std::list<Point>::iterator point = currentPolyline->begin();
//        ++point;
//        while (point != currentPolyline->end()) {
//          snappedOuterRing->addPoint(CGAL::to_double(point->x()), CGAL::to_double(point->y()));
//          ++point;
//        } ++currentPolyline;
//      } snappedPolygon->addRingDirectly(snappedInnerRing);
//    } snappedMultipolygon->OGRGeometryCollection::addGeometryDirectly(snappedPolygon);
//  }
//  
//  return snappedMultipolygon;
//}

//double PolygonRepair::computeRobustness(OGRGeometry *geometry) {
//  double smallestdist = 1e99;
//  
//  if (geometry != NULL) insertConstraints(triangulation, geometry);
//  
//  //-- vertex-vertex distances
//  Vector dist;
//  for (Triangulation::Finite_vertices_iterator curV = triangulation.finite_vertices_begin(); curV != triangulation.finite_vertices_end(); curV++) {
//    // std::cout << "---" << curV->point() << "---" << std::endl;
//    Triangulation::Vertex_circulator vc = triangulation.incident_vertices(curV);
//    Triangulation::Vertex_handle startv = vc;
//    vc++;
//    while (vc != startv) {
//      if (triangulation.is_infinite(vc) == false) {
//        // std::cout << vc->point() << std::endl;
//        dist = curV->point() - vc->point();
//        // std::cout << "distance: " << dist.squared_length()  << std::endl;
//        if (dist.squared_length() < smallestdist)
//          smallestdist = CGAL::to_double(dist.squared_length());
//      }
//      vc++;
//    }
//    dist = curV->point() - startv->point();
//    if (dist.squared_length() < smallestdist)
//      smallestdist = CGAL::to_double(dist.squared_length());
//  }
//  
//  //-- vertex-edge distances
//  for (Triangulation::Finite_vertices_iterator curV = triangulation.finite_vertices_begin(); curV != triangulation.finite_vertices_end(); curV++) {
//    // std::cout << "---" << curV->point() << "---" << std::endl;
//    Triangulation::Face_circulator fc = triangulation.incident_faces(curV);
//    Triangulation::Face_handle startf = fc;
//    fc++;
//    while (fc != startf) {
//      if (triangulation.is_infinite(fc) == false) {
//        int i = fc->index(curV);
//        if (fc->is_constrained(i)) {
//          double d = CGAL::to_double(squared_distance(triangulation.segment(fc, i), fc->vertex(i)->point()));
//          if (d < smallestdist)
//            smallestdist = d;
//        }
//      }
//      fc++;
//    }
//    if (triangulation.is_infinite(startf) == false) {
//      int i = startf->index(curV);
//      if (startf->is_constrained(i)) {
//        double d = CGAL::to_double(squared_distance(triangulation.segment(startf, i), startf->vertex(i)->point()));
//        if (d < smallestdist)
//          smallestdist = d;
//      }
//    }
//  }
//  if (geometry != NULL) triangulation.clear();
//  return smallestdist;
//}

//bool PolygonRepair::saveToShp(OGRGeometry* geometry, const char *fileName) {
//  const char *driverName = "ESRI Shapefile";
//  OGRRegisterAll();
//	OGRSFDriver *driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(driverName);
//	if (driver == NULL) {
//		std::cout << "\tError: OGR Shapefile driver not found." << std::endl;
//		return false;
//	}
//	OGRDataSource *dataSource = driver->Open(fileName, false);
//	if (dataSource != NULL) {
//		std::cout << "\tOverwriting file..." << std::endl;
//		if (driver->DeleteDataSource(dataSource->GetName())!= OGRERR_NONE) {
//			std::cout << "\tError: Couldn't erase file with same name." << std::endl;
//			return false;
//		} OGRDataSource::DestroyDataSource(dataSource);
//	}
//	std::cout << "\tCreating " << fileName << std::endl;
//	dataSource = driver->CreateDataSource(fileName, NULL);
//	if (dataSource == NULL) {
//		std::cout << "\tError: Could not create file." << std::endl;
//		return false;
//	}
//	OGRLayer *layer = dataSource->CreateLayer("polygons", NULL, wkbPolygon, NULL);
//	if (layer == NULL) {
//		std::cout << "\tError: Could not create layer." << std::endl;
//		return false;
//	}
//  OGRFeature *feature = OGRFeature::CreateFeature(layer->GetLayerDefn());
//  // feature->SetField("Name", szName);
//  feature->SetGeometry(geometry);
//  if (layer->CreateFeature(feature) != OGRERR_NONE) {
//    std::cout << "\tError: Could not create feature." << std::endl;
//  }
//  OGRFeature::DestroyFeature(feature);
//  OGRDataSource::DestroyDataSource(dataSource);
//  return true;
//}

void PolygonRepair::insertConstraints(Triangulation &triangulation, MultiPolygon &geometry, bool removeOverlappingConstraints) {
  Triangulation::Vertex_handle va, vb;
  Triangulation::Face_handle faceOfEdge;
  int indexOfEdge;
  
  switch (geometry->getGeometryType()) {
      
    case wkbLineString: {
      OGRLinearRing *ring = static_cast<OGRLinearRing *>(geometry);
      ring->closeRings();
      
      vb = triangulation.insert(Point(ring->getX(0), ring->getY(0)));
      for (int currentPoint = 1; currentPoint < ring->getNumPoints(); ++currentPoint) {
        va = vb;
        vb = triangulation.insert(Point(ring->getX(currentPoint),
                                        ring->getY(currentPoint)),
                                  triangulation.incident_faces(va));
        if (va == vb) continue;
        if (removeOverlappingConstraints && triangulation.is_edge(va, vb, faceOfEdge, indexOfEdge)) {
          if (triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(faceOfEdge, indexOfEdge))) {
            triangulation.insert_constraint(va, vb); // trick to remove a partially overlapping constraint
            triangulation.remove_constraint(va, vb);
          } else triangulation.insert_constraint(va, vb);
        } else triangulation.insert_constraint(va, vb);
      } break;
    }
      
    case wkbPolygon: {
      OGRPolygon *polygon = static_cast<OGRPolygon *>(geometry);
      
      // Outer
      vb = triangulation.insert(Point(polygon->getExteriorRing()->getX(0), polygon->getExteriorRing()->getY(0)));
      for (int currentPoint = 1; currentPoint < polygon->getExteriorRing()->getNumPoints(); ++currentPoint) {
        va = vb;
        vb = triangulation.insert(Point(polygon->getExteriorRing()->getX(currentPoint),
                                        polygon->getExteriorRing()->getY(currentPoint)),
                                  triangulation.incident_faces(va));
        if (va == vb) continue;
        if (removeOverlappingConstraints && triangulation.is_edge(va, vb, faceOfEdge, indexOfEdge)) {
          if (triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(faceOfEdge, indexOfEdge))) {
            triangulation.insert_constraint(va, vb);
            triangulation.remove_constraint(va, vb);
            //std::cout << "Removing constraint <" << va->point() << ", " << vb->point() << ">" << std::endl;
          } else {
            triangulation.insert_constraint(va, vb);
            //std::cout << "Inserting constraint <" << va->point() << ", " << vb->point() << ">" << std::endl;
          }
        } else {
          triangulation.insert_constraint(va, vb);
          //std::cout << "Inserting constraint <" << va->point() << ", " << vb->point() << ">" << std::endl;
        }
      }
      
      // Inner
      for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing) {
        vb = triangulation.insert(Point(polygon->getInteriorRing(currentRing)->getX(0),
                                        polygon->getInteriorRing(currentRing)->getY(0)));
        for (int currentPoint = 1; currentPoint < polygon->getInteriorRing(currentRing)->getNumPoints(); ++currentPoint) {
          va = vb;
          vb = triangulation.insert(Point(polygon->getInteriorRing(currentRing)->getX(currentPoint),
                                          polygon->getInteriorRing(currentRing)->getY(currentPoint)),
                                    triangulation.incident_faces(va));
          if (va == vb) continue;
          if (removeOverlappingConstraints && triangulation.is_edge(va, vb, faceOfEdge, indexOfEdge)) {
            if (triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(faceOfEdge, indexOfEdge))) {
              triangulation.insert_constraint(va, vb);
              triangulation.remove_constraint(va, vb);
              //std::cout << "Removing constraint <" << va->point() << ", " << vb->point() << ">" << std::endl;
            } else {
              triangulation.insert_constraint(va, vb);
              //std::cout << "Inserting constraint <" << va->point() << ", " << vb->point() << ">" << std::endl;
            }
          } else {
            triangulation.insert_constraint(va, vb);
            //std::cout << "Inserting constraint <" << va->point() << ", " << vb->point() << ">" << std::endl;
          }
        }
      } break;
      
    }
      
    case wkbMultiPolygon: {
      OGRMultiPolygon *multipolygon = static_cast<OGRMultiPolygon *>(geometry);
      for (int currentPolygon = 0; currentPolygon < multipolygon->getNumGeometries(); ++currentPolygon) {
        OGRPolygon *polygon = static_cast<OGRPolygon *>(multipolygon->getGeometryRef(currentPolygon));
        
        // Outer
        vb = triangulation.insert(Point(polygon->getExteriorRing()->getX(0),
                                        polygon->getExteriorRing()->getY(0)));
        for (int currentPoint = 1; currentPoint < polygon->getExteriorRing()->getNumPoints(); ++currentPoint) {
          va = vb;
          vb = triangulation.insert(Point(polygon->getExteriorRing()->getX(currentPoint),
                                          polygon->getExteriorRing()->getY(currentPoint)),
                                    triangulation.incident_faces(va));
          if (va == vb) continue;
          if (removeOverlappingConstraints && triangulation.is_edge(va, vb, faceOfEdge, indexOfEdge)) {
            if (triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(faceOfEdge, indexOfEdge))) {
              triangulation.insert_constraint(va, vb);
              triangulation.remove_constraint(va, vb);
            } else triangulation.insert_constraint(va, vb);
          } else triangulation.insert_constraint(va, vb);
        }
        
        // Inner
        for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing) {
          vb = triangulation.insert(Point(polygon->getInteriorRing(currentRing)->getX(0),
                                          polygon->getInteriorRing(currentRing)->getY(0)));
          for (int currentPoint = 1; currentPoint < polygon->getInteriorRing(currentRing)->getNumPoints(); ++currentPoint) {
            va = vb;
            vb = triangulation.insert(Point(polygon->getInteriorRing(currentRing)->getX(currentPoint),
                                            polygon->getInteriorRing(currentRing)->getY(currentPoint)),
                                      triangulation.incident_faces(va));
            if (va == vb) continue;
            if (removeOverlappingConstraints && triangulation.is_edge(va, vb, faceOfEdge, indexOfEdge)) {
              if (triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(faceOfEdge, indexOfEdge))) {
                triangulation.insert_constraint(va, vb);
                triangulation.remove_constraint(va, vb);
              } else triangulation.insert_constraint(va, vb);
            } else triangulation.insert_constraint(va, vb);
          }
        }
      } break;
    }
      
    default:
      std::cerr << "PolygonRepair::insertConstraints: Cannot understand input." << std::endl;
      return;
      break;
  }
  
  // Remove partially even-overlapping subconstraints
  if (!removeOverlappingConstraints) return;
  for (Triangulation::Subconstraint_iterator currentEdge = triangulation.subconstraints_begin();
       currentEdge != triangulation.subconstraints_end();
       ++currentEdge) {
    if (triangulation.number_of_enclosing_constraints(currentEdge->first.first, currentEdge->first.second) % 2 == 0) {
      if (triangulation.is_edge(currentEdge->first.first, currentEdge->first.second, faceOfEdge, indexOfEdge)) {
        if (triangulation.is_constrained(std::pair<Triangulation::Face_handle, int>(faceOfEdge, indexOfEdge))) {
          //std::cout << "Removing constraint <" << currentEdge->first.first->point() << ", " << currentEdge->first.second->point() << ">" << std::endl;
          triangulation.remove_constrained_edge(faceOfEdge, indexOfEdge);
        }
      }
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

//void PolygonRepair::tagPointSet(Triangulation &triangulation, std::list<std::pair<bool, OGRMultiPolygon *> > &geometries) {
//  
//  // Clean tags
//  for (Triangulation::Face_handle currentFace = triangulation.all_faces_begin(); currentFace != triangulation.all_faces_end(); ++currentFace)
//    currentFace->info().clear();
//  
//  std::stack<Triangulation::Face_handle> borderTriangles, taggingStack, untaggingStack;
//  Triangulation::Vertices_in_constraint_iterator currentVertex, nextVertex, lastVertex;
//  Triangulation::Vertex_handle va, vb;
//  Triangulation::Face_handle faceOfSubedge;
//  int indexOfSubedge;
//  bool sameOrder;
//  
//  // Add all repaired outer rings
//  for (std::list<std::pair<bool, OGRMultiPolygon *> >::iterator multipolygon = geometries.begin();
//       multipolygon != geometries.end(); ++multipolygon) {
//    for (int currentPolygon = 0; currentPolygon < multipolygon->second->getNumGeometries(); ++currentPolygon) {
//      OGRPolygon *polygon = static_cast<OGRPolygon *>(multipolygon->second->getGeometryRef(currentPolygon));
//      taggingStack = std::stack<Triangulation::Face_handle>();
//      
//      // Outer
//      if (!multipolygon->first) for (int currentPoint = 0; currentPoint < polygon->getExteriorRing()->getNumPoints(); ++currentPoint) {
//        va = triangulation.insert(Point(polygon->getExteriorRing()->getX(currentPoint),
//                                        polygon->getExteriorRing()->getY(currentPoint)));
//        vb = triangulation.insert(Point(polygon->getExteriorRing()->getX((currentPoint+1)%polygon->getExteriorRing()->getNumPoints()),
//                                        polygon->getExteriorRing()->getY((currentPoint+1)%polygon->getExteriorRing()->getNumPoints())));
//        if (va->point() == vb->point()) continue;
//        currentVertex = triangulation.vertices_in_constraint_begin(va, vb);
//        nextVertex = currentVertex;
//        ++nextVertex;
//        lastVertex = triangulation.vertices_in_constraint_end(va, vb);
//        if (*currentVertex == va) sameOrder = true;
//        else sameOrder = false;
//        while (nextVertex != lastVertex) {
//          if (!sameOrder) {
//            if (!triangulation.is_edge(*currentVertex, *nextVertex, faceOfSubedge, indexOfSubedge)) {
//              std::cerr << "PolygonRepair::tagPointSet: Cannot find edge!" << std::endl;
//              return;
//            }
//          } else {
//            if (!triangulation.is_edge(*nextVertex, *currentVertex, faceOfSubedge, indexOfSubedge)) {
//              std::cerr << "PolygonRepair::tagPointSet: Cannot find edge!" << std::endl;
//              return;
//            }
//          } if (triangulation.number_of_enclosing_constraints(*currentVertex, *nextVertex) % 2 != 0)
//            taggingStack.push(faceOfSubedge);
//          currentVertex = nextVertex;
//          ++nextVertex;
//        }
//      }
//      
//      // Inner
//      else for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing) {
//        for (int currentPoint = 0; currentPoint < polygon->getInteriorRing(currentRing)->getNumPoints(); ++currentPoint) {
//          va = triangulation.insert(Point(polygon->getInteriorRing(currentRing)->getX(currentPoint),
//                                          polygon->getInteriorRing(currentRing)->getY(currentPoint)));
//          vb = triangulation.insert(Point(polygon->getInteriorRing(currentRing)->getX((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints()),
//                                          polygon->getInteriorRing(currentRing)->getY((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints())));
//          if (va->point() == vb->point()) continue;
//          currentVertex = triangulation.vertices_in_constraint_begin(va, vb);
//          nextVertex = currentVertex;
//          ++nextVertex;
//          lastVertex = triangulation.vertices_in_constraint_end(va, vb);
//          if (*currentVertex == va) sameOrder = true;
//          else sameOrder = false;
//          while (nextVertex != lastVertex) {
//            if (sameOrder) {
//              if (!triangulation.is_edge(*currentVertex, *nextVertex, faceOfSubedge, indexOfSubedge)) {
//                std::cerr << "PolygonRepair::tagPointSet: Cannot find edge!" << std::endl;
//                return;
//              }
//            } else {
//              if (!triangulation.is_edge(*nextVertex, *currentVertex, faceOfSubedge, indexOfSubedge)) {
//                std::cerr << "PolygonRepair::tagPointSet: Cannot find edge!" << std::endl;
//                return;
//              }
//            } if (triangulation.number_of_enclosing_constraints(*currentVertex, *nextVertex) % 2 != 0)
//              taggingStack.push(faceOfSubedge);
//            currentVertex = nextVertex;
//            ++nextVertex;
//          }
//        }
//      }
//      
//      // Expand the tags
//      untaggingStack = std::stack<Triangulation::Face_handle>();
//      while (!taggingStack.empty()) {
//        Triangulation::Face_handle currentFace = taggingStack.top();
//        taggingStack.pop();
//        if (currentFace->info().beenTagged()) continue;
//        currentFace->info().isInInterior(true);
//        untaggingStack.push(currentFace);
//        for (int currentEdge = 0; currentEdge < 3; ++currentEdge) {
//          if (!currentFace->neighbor(currentEdge)->info().beenTagged() &&
//              !currentFace->is_constrained(currentEdge)) {
//            taggingStack.push(currentFace->neighbor(currentEdge));
//          }
//        }
//      }
//      
//      // Remove tagged tags
//      while (!untaggingStack.empty()) {
//        untaggingStack.top()->info().beenTagged(false);
//        untaggingStack.pop();
//      }
//    }
//  }
//  
//  // Subtract all repaired inner rings
//  for (std::list<std::pair<bool, OGRMultiPolygon *> >::iterator multipolygon = geometries.begin();
//       multipolygon != geometries.end(); ++multipolygon) {
//    for (int currentPolygon = 0; currentPolygon < multipolygon->second->getNumGeometries(); ++currentPolygon) {
//      OGRPolygon *polygon = static_cast<OGRPolygon *>(multipolygon->second->getGeometryRef(currentPolygon));
//      taggingStack = std::stack<Triangulation::Face_handle>();
//      
//      // Outer
//      if (multipolygon->first) for (int currentPoint = 0; currentPoint < polygon->getExteriorRing()->getNumPoints(); ++currentPoint) {
//        va = triangulation.insert(Point(polygon->getExteriorRing()->getX(currentPoint),
//                                        polygon->getExteriorRing()->getY(currentPoint)));
//        vb = triangulation.insert(Point(polygon->getExteriorRing()->getX((currentPoint+1)%polygon->getExteriorRing()->getNumPoints()),
//                                        polygon->getExteriorRing()->getY((currentPoint+1)%polygon->getExteriorRing()->getNumPoints())));
//        if (va->point() == vb->point()) continue;
//        currentVertex = triangulation.vertices_in_constraint_begin(va, vb);
//        nextVertex = currentVertex;
//        ++nextVertex;
//        lastVertex = triangulation.vertices_in_constraint_end(va, vb);
//        if (*currentVertex == va) sameOrder = true;
//        else sameOrder = false;
//        while (nextVertex != lastVertex) {
//          if (!sameOrder) {
//            if (!triangulation.is_edge(*currentVertex, *nextVertex, faceOfSubedge, indexOfSubedge)) {
//              std::cerr << "PolygonRepair::tagPointSet: Cannot find edge!" << std::endl;
//              return;
//            }
//          } else {
//            if (!triangulation.is_edge(*nextVertex, *currentVertex, faceOfSubedge, indexOfSubedge)) {
//              std::cerr << "PolygonRepair::tagPointSet: Cannot find edge!" << std::endl;
//              return;
//            }
//          } if (triangulation.number_of_enclosing_constraints(*currentVertex, *nextVertex) % 2 != 0) {
//            borderTriangles.push(faceOfSubedge);
//            faceOfSubedge->info().isOnBorder(true);
//            taggingStack.push(faceOfSubedge);
//          } borderTriangles.push(faceOfSubedge->neighbor(indexOfSubedge));
//          faceOfSubedge->neighbor(indexOfSubedge)->info().isOnBorder(true);
//          currentVertex = nextVertex;
//          ++nextVertex;
//        }
//      }
//      
//      // Inner
//      else for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing) {
//        for (int currentPoint = 0; currentPoint < polygon->getInteriorRing(currentRing)->getNumPoints(); ++currentPoint) {
//          va = triangulation.insert(Point(polygon->getInteriorRing(currentRing)->getX(currentPoint),
//                                          polygon->getInteriorRing(currentRing)->getY(currentPoint)));
//          vb = triangulation.insert(Point(polygon->getInteriorRing(currentRing)->getX((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints()),
//                                          polygon->getInteriorRing(currentRing)->getY((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints())));
//          if (va->point() == vb->point()) continue;
//          currentVertex = triangulation.vertices_in_constraint_begin(va, vb);
//          nextVertex = currentVertex;
//          ++nextVertex;
//          lastVertex = triangulation.vertices_in_constraint_end(va, vb);
//          if (*currentVertex == va) sameOrder = true;
//          else sameOrder = false;
//          while (nextVertex != lastVertex) {
//            if (sameOrder) {
//              if (!triangulation.is_edge(*currentVertex, *nextVertex, faceOfSubedge, indexOfSubedge)) {
//                std::cerr << "PolygonRepair::tagPointSet: Cannot find edge!" << std::endl;
//                return;
//              }
//            } else {
//              if (!triangulation.is_edge(*nextVertex, *currentVertex, faceOfSubedge, indexOfSubedge)) {
//                std::cerr << "PolygonRepair::tagPointSet: Cannot find edge!" << std::endl;
//                return;
//              }
//            } if (triangulation.number_of_enclosing_constraints(*currentVertex, *nextVertex) % 2 != 0) {
//              borderTriangles.push(faceOfSubedge);
//              faceOfSubedge->info().isOnBorder(true);
//              taggingStack.push(faceOfSubedge);
//            } borderTriangles.push(faceOfSubedge->neighbor(indexOfSubedge));
//            faceOfSubedge->neighbor(indexOfSubedge)->info().isOnBorder(true);
//            currentVertex = nextVertex;
//            ++nextVertex;
//          }
//        }
//      }
//      
//      // Expand the tags
//      untaggingStack = std::stack<Triangulation::Face_handle>();
//      while (!taggingStack.empty()) {
//        Triangulation::Face_handle currentFace = taggingStack.top();
//        taggingStack.pop();
//        if (currentFace->info().beenTagged()) continue;
//        currentFace->info().isInInterior(false);
//        untaggingStack.push(currentFace);
//        for (int currentEdge = 0; currentEdge < 3; ++currentEdge) {
//          if (!currentFace->neighbor(currentEdge)->info().isOnBorder() &&
//              !currentFace->neighbor(currentEdge)->info().beenTagged()) {
//            taggingStack.push(currentFace->neighbor(currentEdge));
//          }
//        }
//      }
//      
//      // Remove border tags
//      while (!borderTriangles.empty()) {
//        borderTriangles.top()->info().isOnBorder(false);
//        borderTriangles.pop();
//      }
//      
//      // Remove tagged tags
//      while (!untaggingStack.empty()) {
//        untaggingStack.top()->info().beenTagged(false);
//        untaggingStack.pop();
//      }
//    }
//  }
//}

void PolygonRepair::reconstruct(Triangulation &triangulation, MultiPolygon &outPolygons) {
  // std::cout << "Triangulation: " << triangulation.number_of_faces() << " faces, " << triangulation.number_of_vertices() << " vertices." << std::endl;
  if (triangulation.number_of_faces() < 1) {
    return;
  }
  
  // Reconstruct
  for (Triangulation::Finite_faces_iterator seedingFace = triangulation.finite_faces_begin(); seedingFace != triangulation.finite_faces_end(); ++seedingFace) {
    
    if (!seedingFace->info().isInInterior() || seedingFace->info().beenReconstructed()) continue;
    seedingFace->info().beenReconstructed(true);
    if (!seedingFace->info().beenReconstructed()) {
      std::cout << "ERROR! Should be marked as reconstructed!!!" << std::endl;
    }
    
    // Get boundary
    LinearRing vertices;
    if (seedingFace->neighbor(2)->info().isInInterior() && !seedingFace->neighbor(2)->info().beenReconstructed()) {
      seedingFace->neighbor(2)->info().beenReconstructed(true);
      LinearRing l2;
      getBoundary(seedingFace->neighbor(2), seedingFace->neighbor(2)->index(seedingFace), l2);
      vertices.splice(vertices.end(), l2);
    } vertices.push_back(seedingFace->vertex(0));
    if (seedingFace->neighbor(1)->info().isInInterior() && !seedingFace->neighbor(1)->info().beenReconstructed()) {
      seedingFace->neighbor(1)->info().beenReconstructed(true);
      LinearRing l1;
      getBoundary(seedingFace->neighbor(1), seedingFace->neighbor(1)->index(seedingFace), l1);
      vertices.splice(vertices.end(), l1);
    } vertices.push_back(seedingFace->vertex(2));
    if (seedingFace->neighbor(0)->info().isInInterior() && !seedingFace->neighbor(0)->info().beenReconstructed()) {
      seedingFace->neighbor(0)->info().beenReconstructed(true);
      LinearRing l0;
      getBoundary(seedingFace->neighbor(0), seedingFace->neighbor(0)->index(seedingFace), l0);
      vertices.splice(vertices.end(), l0);
    } vertices.push_back(seedingFace->vertex(1));
    
    // Find cutting vertices
    std::set<Triangulation::Vertex_handle> visitedVertices;
    std::set<Triangulation::Vertex_handle> repeatedVertices;
    for (std::list<Triangulation::Vertex_handle>::iterator currentVertex = vertices.begin(); currentVertex != vertices.end(); ++currentVertex) {
      if (!visitedVertices.insert(*currentVertex).second) repeatedVertices.insert(*currentVertex);
    } visitedVertices.clear();
    
    // Cut and join rings in the correct order
    std::list<LinearRing> rings;
    std::list<LinearRing> chainsStack;
    std::set<Triangulation::Vertex_handle> verticesWhereChainsBegin;
    rings.push_back(LinearRing());
    std::list<LinearRing>::reverse_iterator newChain = rings.rbegin();
    for (LinearRing::iterator currentVertex = vertices.begin(); currentVertex != vertices.end(); ++currentVertex) {
      
      // New chain
      if (repeatedVertices.count(*currentVertex) > 0) {
        // Closed by itself
        if (newChain->front() == *currentVertex) {
          // Degenerate (insufficient vertices to be valid)
          if (newChain->size() < 3) {
            std::cerr << "Reconstruction: polygon with less than 3 vertices. Skipped." << std::endl;
            newChain->clear();
          }
          else {
            LinearRing::iterator secondElement = newChain->begin();
            ++secondElement;
            // Degenerate (zero area)
            if (newChain->back() == *secondElement) {
              // std::cerr << "Reconstruction: bridge?. Skipped." << std::endl;
              newChain->clear();
            }
            // Valid
            else {
              rings.push_back(LinearRing());
              newChain = rings.rbegin();
            }
          }
        }
        // Open by itself
        else {
          // Closed with others in stack
          if (verticesWhereChainsBegin.count(*currentVertex)) {
            while (newChain->front() != *currentVertex) {
              newChain->splice(newChain->begin(), chainsStack.back());
              chainsStack.pop_back();
            } verticesWhereChainsBegin.erase(*currentVertex);
            // Degenerate (insufficient vertices to be valid)
            if (newChain->size() < 3) {
              std::cerr << "Reconstruction: polygon with less than 3 vertices. Skipped." << std::endl;
              newChain->clear();
            }
            else {
              std::list<Triangulation::Vertex_handle>::iterator secondElement = newChain->begin();
              ++secondElement;
              // Degenerate (zero area)
              if (newChain->back() == *secondElement) {
                // std::cerr << "Reconstruction: bridge?. Skipped." << std::endl;
                newChain->clear();
              }
              // Valid
              else {
                rings.push_back(LinearRing());
                newChain = rings.rbegin();
              }
            }
          }
          // Open
          else {
            // Not first chain
            if (repeatedVertices.count(newChain->front()) > 0) verticesWhereChainsBegin.insert(newChain->front());
            chainsStack.push_back(LinearRing());
            chainsStack.back().splice(chainsStack.back().begin(), rings.back());
          }
        }
      } newChain->push_back(*currentVertex);
    }
    // Final ring
    while (chainsStack.size() > 0) {
      newChain->splice(newChain->begin(), chainsStack.back());
      chainsStack.pop_back();
    }
    // Degenerate (insufficient vertices to be valid)
    if (newChain->size() < 3) {
      std::cerr << "Reconstruction: polygon with less than 3 vertices. Skipped." << std::endl;
      newChain->clear();
    }
    else {
      std::list<Triangulation::Vertex_handle>::iterator secondElement = newChain->begin();
      ++secondElement;
      // Degenerate (zero area)
      if (newChain->back() == *secondElement) {
        // std::cerr << "Reconstruction: bridge?. Skipped." << std::endl;
        newChain->clear();
      }
      // Valid
      else {
        rings.push_back(LinearRing());
        newChain = rings.rbegin();
      }
    }
    
    // Make rings
    if (rings.size() == 0) continue;
    outPolygons.push_back(Polygon());
    MultiPolygon::reverse_iterator newPolygon = outPolygons.rbegin();
    for (std::list<LinearRing>::iterator currentRing = rings.begin(); currentRing != rings.end(); ++currentRing) {
      if (!currentRing->isClockwise()) {
        newPolygon->push_back(LinearRing());
        Polygon::reverse_iterator newRing = newPolygon->rbegin();
        newRing->splice(newRing->begin(), *currentRing);
        break;
      }
    } for (std::list<LinearRing>::iterator currentRing = rings.begin(); currentRing != rings.end(); ++currentRing) {
      if (currentRing->isClockwise()) {
        newPolygon->push_back(LinearRing());
        Polygon::reverse_iterator newRing = newPolygon->rbegin();
        newRing->splice(newRing->begin(), *currentRing);
      }
    }
  }
}

void PolygonRepair::getBoundary(Triangulation::Face_handle face, int edge, LinearRing &vertices) {
  // Check clockwise edge
  if (face->neighbor(face->cw(edge))->info().isInInterior() && !face->neighbor(face->cw(edge))->info().beenReconstructed()) {
		face->neighbor(face->cw(edge))->info().beenReconstructed(true);
		LinearRing v1;
    getBoundary(face->neighbor(face->cw(edge)), face->neighbor(face->cw(edge))->index(face), v1);
		vertices.splice(vertices.end(), v1);
	}
	
	// Add central vertex
  vertices.push_back(face->vertex(edge));
	
	// Check counterclockwise edge
  if (face->neighbor(face->ccw(edge))->info().isInInterior() && !face->neighbor(face->ccw(edge))->info().beenReconstructed()) {
		face->neighbor(face->ccw(edge))->info().beenReconstructed(true);
		LinearRing v2;
    getBoundary(face->neighbor(face->ccw(edge)), face->neighbor(face->ccw(edge))->index(face), v2);
		vertices.splice(vertices.end(), v2);
	}
}

//void PolygonRepair::printEdges(Triangulation &triangulation) {
//  std::cout << "Edges:" << std::endl;
//  for (Triangulation::Finite_edges_iterator currentEdge = triangulation.finite_edges_begin();
//       currentEdge != triangulation.finite_edges_end();
//       ++currentEdge) {
//    if (triangulation.is_constrained(*currentEdge))
//      std::cout << "\t<" << currentEdge->first->vertex((currentEdge->second+1)%3)->point() << ", " << currentEdge->first->vertex((currentEdge->second+2)%3)->point() << ">" << std::endl;
//  }
//}
//
//void PolygonRepair::printTriangle(Triangulation::Face_handle triangle) {
//  std::cout << "TRIANGLE(" << triangle->vertex(0)->point() << ", " << triangle->vertex(1)->point() << ", " << triangle->vertex(2)->point() << ")";
//}
//
//void PolygonRepair::printChain(LinearRing &chain) {
//  for (LinearRing::iterator currentVertex = chain.begin(); currentVertex != chain.end(); ++currentVertex) {
//    std::cout << (*currentVertex)->point() << ", ";
//  } std::cout << std::endl;
//}
