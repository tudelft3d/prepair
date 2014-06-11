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

#include "Polygon_repair.h"

void usage() {
  std::cout << "=== prepair Help ===\n" << std::endl;
  std::cout << "Usage:   prepair --wkt 'POLYGON((0 0, 0 10, 10 0, 10 10, 0 0))'" << std::endl;
  std::cout << "OR" << std::endl;
  std::cout << "Usage:   prepair -f infile.txt (where infile.txt contains a polygon in WKT)" << std::endl;
  std::cout << "OR" << std::endl;
  std::cout << "Usage:   prepair --ogr infile.shp (first polygon of infile.shp is processed)" << std::endl;
  std::cout << "================================================================================" << std::endl;
  std::cout << "Additional options:" << std::endl;
  std::cout << "--robustness   Compute the robustness of the input/output" << std::endl;
  std::cout << "--setdiff      Uses the point set topology paradigm (default: odd-even paradigm)" << std::endl;
  std::cout << "--min_area AREA Only output polygons larger than AREA" << std::endl;
  std::cout << "--isr GRIDSIZE Snap round the input before repairing" << std::endl;
  std::cout << "--time         Benchmark the different stages of the process" << std::endl;
  
}

int main (int argc, const char * argv[]) {
  
  time_t start_time = time(NULL);
  double min_area = 0;
  double isr_tolerance = 0;
  bool shp_out = false;
  bool compute_robustness = false;
  bool point_set = false;
  bool time_results = false;
  
  if (argc < 2 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
    usage();
    return 0;
  }
  
  OGRGeometry *geometry;
  
  for (int arg_num = 1; arg_num < argc; ++arg_num) {
    //-- time the results
    if (strcmp(argv[arg_num], "--time") == 0) {
      time_results = true;
    }
    
    //-- reading from WKT passed directly
    else if (strcmp(argv[arg_num], "--wkt") == 0) {
      if (arg_num + 1 > argc || argv[arg_num+1][0] == '-') {
        std::cerr << "Error: Invalid parameters" << std::endl;
        return 1;
      } std::string input_wkt(argv[arg_num+1]);
      char *cstr = new char[input_wkt.length()+1];
      std::strcpy(cstr, input_wkt.c_str());
      OGRGeometryFactory::createFromWkt(&cstr, NULL, &geometry);
      if (geometry == NULL) {
        std::cout << "Error: WKT is not valid" << std::endl;
        return 1;
      } ++arg_num;
    }
    
    //-- reading from WKT stored in a text file
    else if (strcmp(argv[arg_num], "-f") == 0) {
      if (arg_num + 1 > argc || argv[arg_num+1][0] == '-') {
        std::cerr << "Error: Invalid parameters" << std::endl;
        return 1;
      } std::ifstream infile;
      infile.open(argv[arg_num+1], std::ios::in | std::ios::binary);
      if (!infile.is_open()) {
        std::cerr << "Error: Could not open file" << std::endl;
        return 1;
      } std::string input_wkt;
      infile.seekg(0, std::ios::end);
      input_wkt.resize(infile.tellg());
      infile.seekg(0, std::ios::beg);
      infile.read(&input_wkt[0], input_wkt.size());
      infile.close();
      char *cstr = new char[input_wkt.length()+1];
      std::strcpy(cstr, input_wkt.c_str());
      OGRGeometryFactory::createFromWkt(&cstr, NULL, &geometry);
      if (geometry == NULL) {
        std::cout << "Error: WKT is not valid" << std::endl;
        return 1;
      } ++arg_num;
    }
    
    //-- reading from a ogr dataset (e.g. shp, geojson, gml, etc.)
    else if (strcmp(argv[arg_num], "--ogr") == 0) {
      OGRRegisterAll();
      if (arg_num + 1 > argc || argv[arg_num+1][0] == '-') {
        std::cerr << "Error: Invalid parameters" << std::endl;
        return 1;
      } OGRDataSource *dataSource = OGRSFDriverRegistrar::Open(argv[arg_num+1], false);
      ++arg_num;
      if (dataSource == NULL) {
        std::cerr << "Error: Could not open file" << std::endl;
        return 1;
      } OGRLayer *dataLayer = dataSource->GetLayer(0);
      dataLayer->ResetReading();
      OGRFeature *feature;
      if (dataLayer->GetFeatureCount(true) > 1) std::cerr << "Warning: More than one feature in input file" << std::endl;
      feature = dataLayer->GetNextFeature();
      if (feature->GetGeometryRef()->getGeometryType() == wkbPolygon)
        geometry = static_cast<OGRPolygon *>(feature->GetGeometryRef());
      else {
        std::cerr << "Error: Feature is not a polygon" << std::endl;
        return 1;
      } OGRFeature::DestroyFeature(feature);
      OGRDataSource::DestroyDataSource(dataSource);
    }
    
    else {
      usage();
      return 0;
    }
  }
  
  if (!time_results) start_time = 0;
  
  Polygon_repair prepair;
  
  char *input_wkt;
  geometry->exportToWkt(&input_wkt);
  std::cout << "Input: " << input_wkt << std::endl;
  
  Multi_polygon<Point> in_polygons, out_polygons;
  Polygon_repair::ogr_to_multi_polygon(geometry, in_polygons);
  prepair.repair_odd_even(in_polygons, out_polygons, start_time);
  
  //-- output well known text
  OGRGeometry *out_geometry = Polygon_repair::multi_polygon_to_ogr(out_polygons);
  char *output_wkt;
  out_geometry->exportToWkt(&output_wkt);
  std::cout << output_wkt << std::endl;
  
  //-- time results
  if (time_results) {
    time_t total_time = time(NULL)-start_time;
    std::cout << "Done! Process finished in " << total_time/60 << " minutes " << total_time%60 << " seconds." << std::endl;
  }
  
  return 0;
}


