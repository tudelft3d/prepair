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

#include "TriangleInfo.h"

TriangleInfo::TriangleInfo() {
  
}

void TriangleInfo::clear() {
  info.reset();
}

bool TriangleInfo::beenTagged() {
  return info[0];
}

void TriangleInfo::beenTagged(bool tagged) {
  info[0] = tagged;
}

bool TriangleInfo::isInInterior() {
  return info[1];
}

void TriangleInfo::isInInterior(bool inInterior) {
  info[0] = true;
  info[1] = inInterior;
}

bool TriangleInfo::isOnBorder() {
  return info[2];
}

void TriangleInfo::isOnBorder(bool onBorder) {
  info[2] = onBorder;
}

bool TriangleInfo::beenReconstructed() {
  return info[3];
}

void TriangleInfo::beenReconstructed(bool reconstructed) {
  info[3] = reconstructed;
}