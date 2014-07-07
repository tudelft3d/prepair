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
#include <boost/program_options.hpp>

int main(int argc, const char *argv[]) {
  
  namespace po = boost::program_options;
  po::options_description main_options("Main options");
  main_options.add_options()
  ("wkt,w", po::value<std::string>()->value_name("'POLYGON(...)'"), "Read WKT passed as a parameter")
  ("wktfile,f", po::value<std::string>()->value_name("PATH"), "Read text file containing one WKT per line")
  ("ogr,i", po::value<std::string>()->value_name("PATH"), "Read file using OGR")
  ("valid,v", "Check if the input is valid")
  ("help,h", "View all options")
  ;
  po::options_description advanced_options("Advanced options");
  advanced_options.add_options()
  ("time,t", "Benchmark the different stages of the process")
  ("setdiff", "Uses the point set paradigm (default: odd-even paradigm)")
  ("minarea", po::value<double>()->value_name("AREA"), "Only output polygons larger than AREA")
  ("isr", po::value<double>()->value_name("GRIDSIZE"), "Snap round the input before repairing")
  ("robustness", "Compute the robustness of the input and output")
  ;
  
  po::options_description all_options;
  all_options.add(main_options);
  all_options.add(advanced_options);
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, all_options), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << "=== prepair help ===" << std::endl;
    std::cout << main_options << std::endl;
    std::cout << advanced_options << std::endl;
    return 0;
  }
  
  time_t start_time = time(NULL);
////  double min_area = 0;
////  bool shp_out = false;
//  bool point_set = false;
//  bool time_results = false;
  
  OGRGeometry *in_geometry;
  OGRDataSource *data_source;
  OGRLayer *data_layer;
  OGRFeature *feature;
  std::ifstream infile;
  
  // Init input
  if (vm.count("wkt")) {
    char *cstr = new char[vm["wkt"].as<std::string>().length()+1];
    std::strcpy(cstr, vm["wkt"].as<std::string>().c_str());
    OGRGeometryFactory::createFromWkt(&cstr, NULL, &in_geometry);
  }
  
  else if (vm.count("wktfile")) {
    infile.open(vm["wktfile"].as<std::string>(), std::ios::in);
    if (!infile.is_open()) {
      std::cerr << "Error: Could not open file" << std::endl;
      return 1;
    }
  }
  
  else if (vm.count("ogr")) {
    OGRRegisterAll();
    data_source = OGRSFDriverRegistrar::Open(vm["ogr"].as<std::string>().c_str(), false);
    if (data_source == NULL) {
      std::cerr << "Error: Could not open file" << std::endl;
      return 1;
    } data_layer = data_source->GetLayer(0);
    data_layer->ResetReading();
  }
  
  else {
    std::cerr << "Error: No input given" << std::endl;
    return 1;
  }
  
  while (true) {

    // Get one polygon
    if (vm.count("wktfile")) {
      if (infile.eof()) {
        infile.close();
        return 0;
      } std::string line;
      std::getline(infile, line);
      char *cstr = new char[line.length()+1];
      OGRGeometryFactory::createFromWkt(&cstr, NULL, &in_geometry);
    }
    
    else if (vm.count("ogr")) {
      feature = data_layer->GetNextFeature();
      if (feature == NULL) {
        OGRDataSource::DestroyDataSource(data_source);
        return 0;
      } in_geometry = feature->GetGeometryRef();
    }
    
    // Do what needs to be done
    Polygon_repair prepair;
    
    if (vm.count("valid")) {
      prepair.is_iso_and_ogc_valid(in_geometry);
    }
    
    OGRGeometry *out_geometry;
    if (vm.count("setdiff")) {
      out_geometry = prepair.repair_point_set(in_geometry);
    } else {
      out_geometry = prepair.repair_odd_even(in_geometry);
    }
    
    // Output results
    char *output_wkt;
    out_geometry->exportToWkt(&output_wkt);
    std::cout << output_wkt << std::endl;
    delete output_wkt;
    
    if (vm.count("wkt")) {
      delete in_geometry;
      return 0;
    }
    
    else if (vm.count("wktfile")) {
      delete in_geometry;
    }
    
    else if (vm.count("ogr")) {
      OGRFeature::DestroyFeature(feature);
    }
  }
  
  // Time results
  if (vm.count("time")) {
    std::time_t total_time = time(NULL)-start_time;
    std::cout << "Done! Process finished in " << total_time/60 << " minutes " << total_time%60 << " seconds." << std::endl;
  }
  
  return 0;
}


