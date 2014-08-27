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

#ifndef EDGE_INFO_H
#define EDGE_INFO_H

#include <CGAL/assertions.h>

class Edge_info {
public:
  Edge_info() {
    CGAL_assertion(sizeof(char) == 1);
    info = 0;
  }
  
  void clear() {
    info = 0;
  }
  
  void add_one() {
    CGAL_precondition(info < 127);
    ++info;
  }
  
  void subtract_one() {
    CGAL_precondition(info > -128);
    --info;
  }
  
  char count() {
    return info;
  }
  
private:
  char info;
};

#endif
