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

#include "Triangle_info.h"

Triangle_info::Triangle_info() {
  info = 0x00;
}

void Triangle_info::clear() {
  info = 0x00;
}

bool Triangle_info::been_tagged() {
  return (info & 0x01) == 0x01;
}

void Triangle_info::been_tagged(bool tagged) {
  if (tagged) info |= 0x01;
  else info &= 0xfe;
}

bool Triangle_info::is_in_interior() {
  return (info & 0x02) == 0x02;
}

void Triangle_info::is_in_interior(bool in_interior) {
  if (in_interior) info |= 0x03;
  else info = (info & 0xfd) | 0x01;
}

bool Triangle_info::is_on_border() {
  return (info & 0x04) == 0x04;
}

void Triangle_info::is_on_border(bool on_border) {
  if (on_border) info |= 0x04;
  else info &= 0xfb;
}

bool Triangle_info::been_reconstructed() {
  return (info & 0x08) == 0x08;
}

void Triangle_info::been_reconstructed(bool reconstructed) {
  if (reconstructed) info |= 0x08;
  else info &= 0xf7;
}