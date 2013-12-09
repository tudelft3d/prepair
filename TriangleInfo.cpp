/*
 Copyright (c) 2009-2013,
 Ken Arroyo Ohori    g.a.k.arroyoohori@tudelft.nl
 Hugo Ledoux         h.ledoux@tudelft.nl
 Martijn Meijers     b.m.meijers@tudelft.nl
 All rights reserved.
 
 // This file is part of prepair: you can redistribute it and/or modify
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
  info = 0x00;
}

void TriangleInfo::clear() {
  info = 0x00;
}

bool TriangleInfo::beenTagged() {
  return (info & 0x01) == 0x01;
}

void TriangleInfo::beenTagged(bool tagged) {
  if (tagged) info |= 0x01;
  else info &= 0xfe;
}

bool TriangleInfo::isInInterior() {
  return (info & 0x02) == 0x02;
}

void TriangleInfo::isInInterior(bool inInterior) {
  if (inInterior) info |= 0x03;
  else info = (info & 0xfd) | 0x01;
}

bool TriangleInfo::isOnBorder() {
  return (info & 0x04) == 0x04;
}

void TriangleInfo::isOnBorder(bool onBorder) {
  if (onBorder) info |= 0x04;
  else info &= 0xfb;
}

bool TriangleInfo::beenReconstructed() {
  return (info & 0x08) == 0x08;
}

void TriangleInfo::beenReconstructed(bool reconstructed) {
  if (reconstructed) info |= 0x08;
  else info &= 0xf7;
}