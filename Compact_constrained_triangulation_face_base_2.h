/*
 Copyright (c) 2014,
 Ken Arroyo Ohori    g.a.k.arroyoohori@tudelft.nl
 All rights reserved.
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 */

#ifndef COMPACT_CONSTRAINED_TRIANGULATION_FACE_BASE_2_H
#define COMPACT_CONSTRAINED_TRIANGULATION_FACE_BASE_2_H

const unsigned char MASKS_P[] = {0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80};
const unsigned char MASKS_N[] = {0xfe, 0xfd, 0xfb, 0xf7, 0xef, 0xdf, 0xbf, 0x7f};

#include <CGAL/Triangulation_face_base_2.h>

template <class GT, class FB = CGAL::Triangulation_face_base_2<GT> >
class Compact_constrained_triangulation_face_base_2 : public FB {
public:
  typedef typename FB::Triangulation_data_structure Triangulation_data_structure;
  typedef typename Triangulation_data_structure::Vertex_handle Vertex_handle;
  typedef typename Triangulation_data_structure::Face_handle Face_handle;
  
  template <class TDS_2>
  struct Rebind_TDS {
    typedef Compact_constrained_triangulation_face_base_2<GT, typename FB::template Rebind_TDS<TDS_2>::Other> Other;
  };
protected:
  unsigned char constrained;
  
public:
  Compact_constrained_triangulation_face_base_2() : FB() {
    constrained = 0x00;
  }
  
  Compact_constrained_triangulation_face_base_2(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2) : FB(v0, v1, v2) {
    constrained = 0x00;
  }
  
  Compact_constrained_triangulation_face_base_2(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Face_handle n0, Face_handle n1, Face_handle n2) : FB(v0, v1, v2, n0, n1, n2) {
    constrained = 0x00;
  }
  
  bool is_constrained(int i) const {
    CGAL_triangulation_precondition( i == 0 || i == 1 || i == 2);
    return (constrained & MASKS_P[i]) == MASKS_P[i];
  }
  
  void set_constraint(int i, bool b) {
    CGAL_triangulation_precondition( i == 0 || i == 1 || i == 2);
    if (b) constrained |= MASKS_P[i];
    else constrained &= MASKS_N[i];
  }
  
  void set_constraints(bool c0, bool c1, bool c2) {
    set_constraint(0, c0);
    set_constraint(1, c1);
    set_constraint(2, c2);
  }
  
  void reorient() {
    FB::reorient();
    set_constraints(is_constrained(1), is_constrained(0), is_constrained(2));
  }
  
  void ccw_permute() {
    FB::ccw_permute();
    set_constraints(is_constrained(2), is_constrained(0), is_constrained(1));
  }
  
  void cw_permute() {
    FB::cw_permute();
    set_constraints(is_constrained(1), is_constrained(2), is_constrained(0));
  }
};

#endif