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

#ifndef TRIANGLEINFO_H
#define TRIANGLEINFO_H

#include <bitset>

class Triangle_info {
public:
  Triangle_info();
  
  void clear();
  bool been_tagged();
  void been_tagged(bool tagged);
  bool is_in_interior();
  void is_in_interior(bool in_interior);
  bool is_on_border();
  void is_on_border(bool on_order);
  bool been_reconstructed();
  void been_reconstructed(bool reconstructed);
  
private:
  std::bitset<4> info;
};

#endif
