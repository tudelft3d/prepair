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

class Edge_info {
public:
  Edge_info() {
    info = 0x00;
  }
  
  void clear() {
    info = 0x00;
  }
  
  bool to_fill_in() {
    return (info & 0x01) == 0x01;
  }
  
  void to_fill_in(bool fill_in) {
    if (fill_in) info |= 0x01;
    else info &= 0xfe;
  }
  
  bool to_carve_out() {
    return (info & 0x02) == 0x02;
  }
  
  void to_carve_out(bool carve_out) {
    if (carve_out) info |= 0x02;
    else info &= 0xfd;
  }
  
private:
  unsigned char info;
};

#endif
