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

#include <list>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_repair/repair.h>
#include <CGAL/linear_least_squares_fitting_3.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Point_3 = Kernel::Point_3;
using Plane_3 = Kernel::Plane_3;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel>;

struct Polygon_repair {
  OGRGeometry *geometry;

  void repair() {
    if (!geometry || geometry->IsEmpty()) {
      geometry = new OGRPolygon();
      return;
    }

    const bool is_3d = geometry->Is3D();
    if (is_3d) compute_plane();

    Multipolygon_with_holes_2 cgal_input = to_cgal_multipolygon(geometry, is_3d);

    Multipolygon_with_holes_2 cgal_output;
    if (cgal_input.number_of_polygons_with_holes() > 0)
      cgal_output = CGAL::Polygon_repair::repair(cgal_input);

    OGRGeometry *new_geom = from_cgal_multipolygon(cgal_output, is_3d);
    delete geometry;
    geometry = new_geom;
  }

private:
  Plane_3 best_plane;

  void compute_plane() {
    std::list<Point_3> points;
    extract_points(geometry, points);
    CGAL::linear_least_squares_fitting_3(
      points.begin(), points.end(), best_plane, CGAL::Dimension_tag<0>());
  }

  static void extract_points(OGRGeometry *g, std::list<Point_3> &points) {
    switch (g->getGeometryType()) {
      case wkbLineString25D: {
        auto *ring = static_cast<OGRLinearRing *>(g);
        ring->closeRings();
        for (int i = 1; i < ring->getNumPoints(); ++i)
          points.emplace_back(ring->getX(i), ring->getY(i), ring->getZ(i));
        break;
      }
      case wkbPolygon25D: {
        auto *polygon = static_cast<OGRPolygon *>(g);
        extract_points(polygon->getExteriorRing(), points);
        for (int i = 0; i < polygon->getNumInteriorRings(); ++i)
          extract_points(polygon->getInteriorRing(i), points);
        break;
      }
      case wkbMultiPolygon25D: {
        auto *mp = static_cast<OGRMultiPolygon *>(g);
        for (int i = 0; i < mp->getNumGeometries(); ++i)
          extract_points(mp->getGeometryRef(i), points);
        break;
      }
      default: break;
    }
  }

  static Point_2 point_to_2d(OGRLinearRing *ring, int i) {
    return Point_2(ring->getX(i), ring->getY(i));
  }

  Point_2 point_to_2d_projected(OGRLinearRing *ring, int i) const {
    return best_plane.to_2d(Point_3(ring->getX(i), ring->getY(i), ring->getZ(i)));
  }

  Polygon_2 ring_to_cgal(OGRLinearRing *ring, bool is_3d) {
    ring->closeRings();
    Polygon_2 poly;
    poly.reserve(ring->getNumPoints());
    for (int i = 0; i < ring->getNumPoints(); ++i)
      poly.push_back(is_3d ? point_to_2d_projected(ring, i) : point_to_2d(ring, i));
    return poly;
  }

  Polygon_with_holes_2 polygon_to_cgal(OGRPolygon *polygon, bool is_3d) {
    Polygon_2 outer = ring_to_cgal(polygon->getExteriorRing(), is_3d);
    std::vector<Polygon_2> holes;
    holes.reserve(polygon->getNumInteriorRings());
    for (int i = 0; i < polygon->getNumInteriorRings(); ++i)
      holes.push_back(ring_to_cgal(polygon->getInteriorRing(i), is_3d));
    return Polygon_with_holes_2(std::move(outer), holes.begin(), holes.end());
  }

  Multipolygon_with_holes_2 to_cgal_multipolygon(OGRGeometry *g, bool is_3d) {
    Multipolygon_with_holes_2 mp;
    switch (g->getGeometryType()) {
      case wkbPolygon:
      case wkbPolygon25D:
        mp.add_polygon_with_holes(polygon_to_cgal(static_cast<OGRPolygon *>(g), is_3d));
        break;
      case wkbMultiPolygon:
      case wkbMultiPolygon25D: {
        auto *multi = static_cast<OGRMultiPolygon *>(g);
        for (int i = 0; i < multi->getNumGeometries(); ++i)
          mp.add_polygon_with_holes(
            polygon_to_cgal(static_cast<OGRPolygon *>(multi->getGeometryRef(i)), is_3d));
        break;
      }
      default: break;
    }
    return mp;
  }

  OGRLinearRing *cgal_to_ring(const Polygon_2 &poly, bool is_3d) const {
    auto *ring = new OGRLinearRing();
    for (auto it = poly.vertices_begin(); it != poly.vertices_end(); ++it) {
      if (is_3d) {
        Point_3 p3 = best_plane.to_3d(*it);
        ring->addPoint(p3.x(), p3.y(), p3.z());
      } else {
        ring->addPoint(it->x(), it->y());
      }
    }
    ring->closeRings();
    return ring;
  }

  OGRGeometry *from_cgal_multipolygon(const Multipolygon_with_holes_2 &mp, bool is_3d) {
    if (mp.number_of_polygons_with_holes() == 0)
      return new OGRPolygon();

    auto *result = new OGRMultiPolygon();
    for (auto it = mp.polygons_with_holes().begin(); it != mp.polygons_with_holes().end(); ++it) {
      auto *poly = new OGRPolygon();
      poly->addRingDirectly(cgal_to_ring(it->outer_boundary(), is_3d));
      for (auto hit = it->holes_begin(); hit != it->holes_end(); ++hit)
        poly->addRingDirectly(cgal_to_ring(*hit, is_3d));
      result->addGeometryDirectly(poly);
    }

    if (result->getNumGeometries() == 1) {
      OGRPolygon *single = static_cast<OGRPolygon *>(result->getGeometryRef(0)->clone());
      delete result;
      return single;
    }
    return result;
  }
};

#endif
