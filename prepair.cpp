// prepair
//
// Copyright © 2009-2020,
// Ken Arroyo Ohori    k.ohori@tudelft.nl
// Hugo Ledoux         h.ledoux@tudelft.nl
// Martijn Meijers     b.m.meijers@tudelft.nl
// All rights reserved.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <ogrsf_frmts.h>
#include "Polygon_repair.h"

int main(int argc, const char *argv[]) {
  
  boost::program_options::options_description input_options("input options");
  input_options.add_options()
    ("wkt,w", boost::program_options::value<std::string>()->value_name("'POLYGON(...)'"), "Read WKT passed directly as a parameter")
    ("wktfile,t", boost::program_options::value<std::string>()->value_name("PATH"), "Read text file containing one WKT per line")
    ("ogrin,i", boost::program_options::value<std::string>()->value_name("path"), "Read data source using OGR")
    ("help,h", "View all options");
  boost::program_options::options_description output_options("output options");
  output_options.add_options()
    ("ogrout,o", boost::program_options::value<std::string>()->value_name("path"), "Output to a file using OGR");
  boost::program_options::options_description all_options;
  all_options.add(input_options).add(output_options);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::parse_command_line(argc, argv, all_options), vm);
  boost::program_options::notify(vm);
  
  GDALDataset *in_dataset = NULL;
  std::ifstream in_file;
  OGRGeometry *in_geometry = NULL;
  
  // Help or invalid args
  if (argc < 2 || vm.count("help")) {
    std::cout << "usage: prepair input_option [output_option]" << std::endl;
    std::cout << "if no output option is provided, prepair will output WKT directly" << std::endl << std::endl;
    std::cout << input_options << std::endl;
    std::cout << output_options << std::endl;
    return 0;
  }
  
  // Input
  if (vm.count("wkt")) {
    std::string wkt = vm["wkt"].as<std::string>();
    OGRErr error = OGRGeometryFactory::createFromWkt(wkt.c_str(), NULL, &in_geometry);
    if (error != OGRERR_NONE) {
      switch (error) {
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
      } return 1;
    } if (in_geometry->IsEmpty()) {
      std::cerr << "Error: empty geometry" << std::endl;
      return 1;
    } if ((in_geometry->getGeometryType() != wkbPolygon) &&
          (in_geometry->getGeometryType() != wkbPolygon25D) &&
          (in_geometry->getGeometryType() != wkbMultiPolygon) &&
          (in_geometry->getGeometryType() != wkbMultiPolygon25D)) {
       std::cerr << "Error: geometry must be Polygon or MultiPolygon (in 2D or 3D)" << std::endl;
       return 1;
     }
  }
  
  else if (vm.count("wktfile")) {
    in_file.open(vm["wktfile"].as<std::string>().c_str(), std::ios::in);
    if (!in_file.is_open()) {
      std::cerr << "Error: couldn't open file" << std::endl;
      return 1;
    }
  }
  
  else if (vm.count("ogrin")) {
    GDALAllRegister();
    in_dataset = (GDALDataset *)GDALOpenEx(vm["ogrin"].as<std::string>().c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL);
    if (in_dataset == NULL) {
      std::cerr << "Error: couldn't open file" << std::endl;
      return 1;
    }
  }
  
  // Output
  if (vm.count("ogrout")) {
    GDALAllRegister();
    boost::filesystem::path extension = boost::filesystem::path(vm["ogrout"].as<std::string>()).extension();
    std::string out_driver_name;
    if (extension.compare(".csv") == 0) out_driver_name = "CSV";
    else if (extension.compare(".dxf") == 0) out_driver_name = "DXF";
    else if (extension.compare(".gdb") == 0) out_driver_name = "FileGDB";
    else if (extension.compare(".json") == 0) out_driver_name = "GeoJSON";
    else if (extension.compare(".geojson") == 0) out_driver_name = "GeoJSON";
    else if (extension.compare(".gml") == 0) out_driver_name = "GML";
    else if (extension.compare(".gpkg") == 0) out_driver_name = "GPKG";
    else if (extension.compare(".kml") == 0) out_driver_name = "KML";
    else if (extension.compare(".shp") == 0) out_driver_name = "ESRI Shapefile";
    else {
      std::cout << "Error: unknown output format" << std::endl;
      return 1;
    } GDALDriver *out_driver = GetGDALDriverManager()->GetDriverByName(out_driver_name.c_str());
    if (out_driver == NULL) {
      std::cerr << "Error: OGR driver not found" << std::endl;
      return 1;
    }
    
    if (boost::filesystem::exists(boost::filesystem::path(vm["ogrout"].as<std::string>()))) {
      std::cout << "Overwriting " << out_driver_name << " file " << vm["ogrout"].as<std::string>() << "..." << std::endl;
      if (out_driver->Delete(vm["ogrout"].as<std::string>().c_str())!= OGRERR_NONE) {
        std::cerr << "Error: couldn't overwrite file" << std::endl;
        return 1;
      }
    } else std::cout << "Writing " << out_driver_name << " file " << vm["ogrout"].as<std::string>() << "..." << std::endl;
    GDALDataset *out_dataset = out_driver->Create(vm["ogrout"].as<std::string>().c_str(), 0, 0, 0, GDT_Unknown, NULL);
    if (out_dataset == NULL) {
      std::cout << "Error: couldn't create file" << std::endl;
      return 1;
    }
    
    if (vm.count("wkt")) {
      OGRLayer *out_layer = out_dataset->CreateLayer("polygons");
      if (out_layer == NULL) {
        std::cout << "Error: couldn't create layer." << std::endl;
        return 1;
      } OGRFeature *out_feature = OGRFeature::CreateFeature(out_layer->GetLayerDefn());
      Polygon_repair pr;
      pr.geometry = in_geometry;
      pr.repair();
      OGRGeometry *out_geometry = pr.geometry;
      out_feature->SetGeometry(out_geometry);
      if (out_layer->CreateFeature(out_feature) != OGRERR_NONE) {
        std::cout << "Error: couldn't create feature." << std::endl;
        return 1;
      }
    }
    
    else if (vm.count("wktfile")) {
      OGRLayer *out_layer = out_dataset->CreateLayer("polygons");
      if (out_layer == NULL) {
        std::cout << "Error: couldn't create layer." << std::endl;
        return 1;
      } std::string line;
      while (std::getline(in_file, line)) {
        OGRErr error = OGRGeometryFactory::createFromWkt(line.c_str(), NULL, &in_geometry);
        if (error != OGRERR_NONE) {
          OGRFeature *out_feature = OGRFeature::CreateFeature(out_layer->GetLayerDefn());
          Polygon_repair pr;
          pr.geometry = in_geometry;
          pr.repair();
          OGRGeometry *out_geometry = pr.geometry;
          out_feature->SetGeometry(out_geometry);
          if (out_layer->CreateFeature(out_feature) != OGRERR_NONE) {
            std::cout << "Error: couldn't create feature." << std::endl;
            return 1;
          }
        }
      }
    }
    
    else if (vm.count("ogrin")) {
      for (OGRLayer *layer: in_dataset->GetLayers()) {
        OGRLayer *out_layer = out_dataset->CreateLayer(layer->GetName());
        if (out_layer == NULL) {
          std::cout << "Error: couldn't create layer." << std::endl;
          return 1;
        } for (OGRFeatureUniquePtr &feature: *layer) {
          OGRFeature *out_feature = feature->Clone();
          Polygon_repair pr;
          pr.geometry = feature->GetGeometryRef();
          pr.repair();
          OGRGeometry *out_geometry = pr.geometry;
          out_feature->SetGeometry(out_geometry);
          if (out_layer->CreateFeature(out_feature) != OGRERR_NONE) {
            std::cout << "Error: couldn't create feature." << std::endl;
            return 1;
          }
        }
      }
    }
    
    GDALClose(out_dataset);
  }
  
  else {
    
    if (vm.count("wkt")) {
      Polygon_repair pr;
      pr.geometry = in_geometry;
      pr.repair();
      OGRGeometry *out_geometry = pr.geometry;
      char *output_wkt;
      out_geometry->exportToWkt(&output_wkt);
      std::cout << output_wkt << std::endl;
      delete output_wkt;
    }
    
    else if (vm.count("wktfile")) {
      std::string line;
      while (std::getline(in_file, line)) {
        OGRErr error = OGRGeometryFactory::createFromWkt(line.c_str(), NULL, &in_geometry);
        if (error != OGRERR_NONE) {
          Polygon_repair pr;
          pr.geometry = in_geometry;
          pr.repair();
          OGRGeometry *out_geometry = pr.geometry;
          char *output_wkt;
          out_geometry->exportToWkt(&output_wkt);
          std::cout << output_wkt << std::endl;
          delete output_wkt;
        } else std::cout << std::endl;
      }
    }
    
    else if (vm.count("ogrin")) {
      for (OGRLayer *layer: in_dataset->GetLayers()) {
        for (OGRFeatureUniquePtr &feature: *layer) {
          Polygon_repair pr;
          pr.geometry = feature->GetGeometryRef();
          pr.repair();
          OGRGeometry *out_geometry = pr.geometry;
          char *output_wkt;
          out_geometry->exportToWkt(&output_wkt);
          std::cout << output_wkt << std::endl;
          delete output_wkt;
        }
      }
    }
    
  }
  
  return 0;
}
