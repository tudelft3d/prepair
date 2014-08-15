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

#ifndef TRIANGULATION_FACE_BASE_WITH_INFO_ON_FACE_AND_HALFEDGES_2_H
#define TRIANGULATION_FACE_BASE_WITH_INFO_ON_FACE_AND_HALFEDGES_2_H

template <class FI, class HEI, class GT, class FB = CGAL::Triangulation_face_base_2<GT> >
class Triangulation_face_base_with_info_on_face_and_halfedges_2 : public FB {
public:
  typedef typename FB::Triangulation_data_structure Triangulation_data_structure;
  typedef typename Triangulation_data_structure::Vertex_handle Vertex_handle;
  typedef typename Triangulation_data_structure::Face_handle Face_handle;
  typedef FI Face_info;
  typedef HEI Halfedge_info;
  
  template <class TDS_2>
  struct Rebind_TDS {
    typedef Triangulation_face_base_with_info_on_face_and_halfedges_2<FI, HEI, GT, typename FB::template Rebind_TDS<TDS_2>::Other> Other;
  };
  
protected:
  Face_info fi;
  Halfedge_info ei[3];
  
public:
  Triangulation_face_base_with_info_on_face_and_halfedges_2() : FB() {}
  Triangulation_face_base_with_info_on_face_and_halfedges_2(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2) : FB(v0, v1, v2) {}
  Triangulation_face_base_with_info_on_face_and_halfedges_2(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Face_handle n0, Face_handle n1, Face_handle n2) : FB(v0, v1, v2, n0, n1, n2) {}
  
  const Face_info &info() const {
    return fi;
  }
  
  Face_info &info() {
    return fi;
  }
  
  const Halfedge_info &halfedge_info(int i) const {
    return ei[i];
  }
  
  Halfedge_info &halfedge_info(int i) {
    return ei[i];
  }
};

#endif