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

// TODO: Maybe create a marker manager in the future
class TriangleInfo {
public:
  TriangleInfo();
  
  void clear();
  bool beenTagged();
  void beenTagged(bool tagged);
  bool isInInterior();
  void isInInterior(bool inInterior);
  bool isOnBorder();
  void isOnBorder(bool onBorder);
  bool beenReconstructed();
  void beenReconstructed(bool reconstructed);
  
private:
  std::bitset<4> info;
};

#endif
