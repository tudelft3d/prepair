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
  ("wkt,w", po::value<std::string>()->value_name("'POLYGON(...)'"), "Read WKT passed directly as a parameter")
  ("wktfile,f", po::value<std::string>()->value_name("PATH"), "Read text file containing one WKT per line")
  ("ogr,i", po::value<std::string>()->value_name("PATH"), "Read file using OGR")
  ("help,h", "View all options")
  ;
  po::options_description advanced_options("Advanced options");
  advanced_options.add_options()
  ("setdiff", "Uses the point set paradigm (default: odd-even paradigm)")
  ("minarea", po::value<double>()->value_name("AREA"), "Only output polygons larger than AREA")
  ("shpOut", po::value<std::string>()->value_name("PATH"), "Output to a shapefile")
  ("time,t", "Benchmark the different stages of the process")
  ;
  po::options_description hidden_options("Hidden options");
  hidden_options.add_options()
  ("valid,v", "Check if the input is valid rather than repair")
  ("robustness", "Compute the robustness of the input and output")
  ("noOut", "Compute but do not print the output")
  ;
  
  po::options_description all_options;
  all_options.add(main_options).add(advanced_options).add(hidden_options);
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, all_options), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << "=== prepair help ===" << std::endl;
    std::cout << main_options << std::endl;
    std::cout << advanced_options << std::endl;
    return 0;
  }
  
  std::clock_t start_time = clock();
  bool time_results = false;
  
  OGRGeometry *in_geometry;
  OGRDataSource *in_data_source;
  OGRLayer *in_data_layer;
  OGRFeature *in_feature;
  std::ifstream in_file;
  
  // Init input
  if (vm.count("wkt")) {
    char *cstr = new char[vm["wkt"].as<std::string>().length()+1];
    std::strcpy(cstr, vm["wkt"].as<std::string>().c_str());
    OGRErr err = OGRGeometryFactory::createFromWkt(&cstr, NULL, &in_geometry);
    if (err != OGRERR_NONE) {
      switch (err) {
        case OGRERR_UNSUPPORTED_GEOMETRY_TYPE:
          std::cerr << "Error: unsupported geometry" << std::endl;
          break;
        case OGRERR_NOT_ENOUGH_DATA:
          std::cerr << "Error: degenerate input" << std::endl;
          break;
        case OGRERR_CORRUPT_DATA:
          std::cerr << "Error: corrupted input" << std::endl;
          break;
        default:
          std::cerr << "Error" << std::endl;
          break;
      }
      return 1;
    }
    if (in_geometry->IsEmpty() == 1) {
      std::cerr << "Error: empty geometry" << std::endl;
      return 1;
    }
    if ( (in_geometry->getGeometryType() != wkbPolygon) && 
         (in_geometry->getGeometryType() != wkbMultiPolygon) ) {
      std::cerr << "Error: geometry must be Polygon or MultiPolygon" << std::endl;
      return 1;
    }
  }
  
  else if (vm.count("wktfile")) {
    in_file.open(vm["wktfile"].as<std::string>().c_str(), std::ios::in);
    if (!in_file.is_open()) {
      std::cerr << "Error: Could not open file" << std::endl;
      return 1;
    } else {
      std::cout << "Opened: " << vm["wktfile"].as<std::string>() << std::endl;
    }
  }
  
  else if (vm.count("ogr")) {
    OGRRegisterAll();
    in_data_source = OGRSFDriverRegistrar::Open(vm["ogr"].as<std::string>().c_str(), false);
    if (in_data_source == NULL) {
      std::cerr << "Error: Could not open file" << std::endl;
      return 1;
    } in_data_layer = in_data_source->GetLayer(0);
    in_data_layer->ResetReading();
  }
  
  else {
    std::cerr << "Error: No input given" << std::endl;
    std::cout << "=== prepair help ===" << std::endl;
    std::cout << main_options << std::endl;
    std::cout << advanced_options << std::endl;
    return 1;
  }
  
  if (vm.count("time")) {
    time_results = true;
  }
  
  OGRSFDriver *out_driver;
  OGRDataSource *out_data_source;
  OGRLayer *out_layer;
  if (vm.count("shpOut")) {
    const char *out_driver_name = "ESRI Shapefile";
    OGRRegisterAll();
    out_driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(out_driver_name);
    if (out_driver == NULL) {
      std::cout << "Error: OGR Shapefile driver not found." << std::endl;
      return false;
    } out_data_source = out_driver->Open(vm["shpOut"].as<std::string>().c_str(), false);
    if (out_data_source != NULL) {
      std::cout << "Overwriting " << vm["shpOut"].as<std::string>() << "..." << std::endl;
      if (out_driver->DeleteDataSource(out_data_source->GetName())!= OGRERR_NONE) {
        std::cout << "Error: Couldn't erase file with same name." << std::endl;
        return false;
      } OGRDataSource::DestroyDataSource(out_data_source);
    } else {
      std::cout << "Creating " << vm["shpOut"].as<std::string>() << "..." << std::endl;
    } out_data_source = out_driver->CreateDataSource(vm["shpOut"].as<std::string>().c_str(), NULL);
    if (out_data_source == NULL) {
      std::cout << "Error: Could not create file." << std::endl;
      return false;
    } out_layer = out_data_source->CreateLayer("polygons", NULL, wkbPolygon, NULL);
    if (out_layer == NULL) {
      std::cout << "Error: Could not create layer." << std::endl;
      return false;
    }
  }
  
  while (true) {

    // Get one polygon (WKT file)
    if (vm.count("wktfile")) {
      if (in_file.eof()) {
        in_file.close();
        break;
      } std::string line;
      std::getline(in_file, line);
      char *cstr = new char[line.length()+1];
      std::strcpy(cstr, line.c_str());
      OGRGeometryFactory::createFromWkt(&cstr, NULL, &in_geometry);
    }
    
    // Get one polygon (OGR)
    else if (vm.count("ogr")) {
      in_feature = in_data_layer->GetNextFeature();
      if (in_feature == NULL) {
        OGRDataSource::DestroyDataSource(in_data_source);
        break;
      } in_geometry = in_feature->GetGeometryRef();
    }
    
    if (in_geometry == NULL) break;
    
    // Do what needs to be done
    Polygon_repair prepair;
    
    if (vm.count("valid")) {
      prepair.is_iso_and_ogc_valid(in_geometry);
    }
    
    OGRGeometry *out_geometry;
    if (vm.count("setdiff")) {
      out_geometry = prepair.repair_point_set(in_geometry, time_results);
    } else {
      out_geometry = prepair.repair_odd_even(in_geometry, time_results);
    }
    
    // Remove small parts
    if (vm.count("minarea")) {
      prepair.remove_small_parts(out_geometry, vm["minarea"].as<double>());
    }
    
    // Output results
    if (vm.count("shpOut")) {
      OGRFeature *out_feature = OGRFeature::CreateFeature(out_layer->GetLayerDefn());
      out_feature->SetGeometry(out_geometry);
      if (out_layer->CreateFeature(out_feature) != OGRERR_NONE) {
        std::cout << "Error: Could not create feature." << std::endl;
      } OGRFeature::DestroyFeature(out_feature);
    } else if (!vm.count("noOut")) {
      char *output_wkt;
      out_geometry->exportToWkt(&output_wkt);
      std::cout << output_wkt << std::endl;
      delete output_wkt;
    }
    
    if (vm.count("wkt")) {
      delete in_geometry;
      break;
    }
    
    else if (vm.count("wktfile")) {
      delete in_geometry;
    }
    
    else if (vm.count("ogr")) {
      OGRFeature::DestroyFeature(in_feature);
    }
  }
  
  if (vm.count("shpOut")) {
    OGRDataSource::DestroyDataSource(out_data_source);
  }
  
  // Time results
  if (time_results) {
    std::clock_t total_time = clock()-start_time;
    std::cout << "Done! Process finished in " << total_time / double(CLOCKS_PER_SEC) << " seconds." << std::endl;
  }
  
  return 0;
}


