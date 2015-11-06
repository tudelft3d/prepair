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

#include "PolygonRepair.h"

double  minArea = 0;
double  isrTolerance = 0;
bool    shpOut = false;
bool    computeRobustness = false;
bool    pointSet = false;
bool    timeResults = false;

void    usage();
bool    savetoshp(OGRMultiPolygon* multiPolygon);


int main (int argc, const char * argv[]) {
  
  time_t startTime = time(NULL);
  
  if (argc < 2 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
    usage();
    return(0);
  }
  
  OGRGeometry *geometry;
  
  for (int argNum = 1; argNum < argc; ++argNum) {
    //-- whether to compute the robustness or not
    if (strcmp(argv[argNum], "--robustness") == 0) {
      computeRobustness = true;
    }
    
    //-- whether to use the point set topology paradigm or not
    //-- if not, odd-even paradigm is used by default
    else if (strcmp(argv[argNum], "--setdiff") == 0) {
      pointSet = true;
    }
    
    //-- mininum area to keep in output
    else if (strcmp(argv[argNum], "--minarea") == 0) {
      minArea = atof(argv[argNum+1]);
      ++argNum;
    }
    
    //-- ISR snapping tolerance
    else if (strcmp(argv[argNum], "--isr") == 0) {
      isrTolerance = atof(argv[argNum+1]);
      ++argNum;
      // TODO : scale dataset if tolerance < 1.0 because of CGAL bug
    }
    
    //-- output a shapefile (out.shp) instead of a WKT
    else if (strcmp(argv[argNum], "--shpOut") == 0) {
      shpOut = true;
    }
    
    //-- time the results
    else if (strcmp(argv[argNum], "--time") == 0) {
      timeResults = true;
    }
    
    //-- reading from WKT passed directly
    else if (strcmp(argv[argNum], "--wkt") == 0) {
      unsigned int bufferSize = 100000000;
      char *inputWKT = (char *)malloc(bufferSize*sizeof(char));
      strcpy(inputWKT, argv[argNum+1]);
      ++argNum;
      OGRErr err = OGRGeometryFactory::createFromWkt(&inputWKT, NULL, &geometry);
      if (err != OGRERR_NONE) {
        switch (err) {
          case OGRERR_UNSUPPORTED_GEOMETRY_TYPE:
            std::cerr << "Error: geometry must be Polygon or MultiPolygon" << std::endl;
            break;
          case OGRERR_NOT_ENOUGH_DATA:
          case OGRERR_CORRUPT_DATA:
            std::cerr << "Error: corrupted input" << std::endl;
            break;
          default:
            std::cerr << "Error: corrupted input" << std::endl;
            break;
        }
        return 1;
      }
      if (geometry->IsEmpty() == 1) {
        std::cerr << "Error: empty geometry" << std::endl;
        return 1;
      }
      if ( (geometry->getGeometryType() != wkbPolygon) && 
           (geometry->getGeometryType() != wkbMultiPolygon) ) {
        std::cerr << "Error: geometry must be Polygon or MultiPolygon" << std::endl;
        return 1;
      }
    }
    
    //-- reading from WKT stored in first line of a text file
    else if (strcmp(argv[argNum], "-f") == 0) {
      unsigned int bufferSize = 100000000;
      char *inputWKT = (char *)malloc(bufferSize*sizeof(char));
      if (argNum + 1 <= argc - 1 && argv[argNum+1][0] != '-') {
        std::ifstream infile(argv[argNum+1], std::ifstream::in);
        infile.getline(inputWKT, bufferSize);
        ++argNum;
      } else {
        std::cerr << "Error: Missing input file name." << std::endl;
        return 1;
      }
      OGRGeometryFactory::createFromWkt(&inputWKT, NULL, &geometry);
      if (geometry == NULL) {
        std::cout << "Error: WKT is not valid" << std::endl;
        return 1;
      }
    }
    
    //-- reading from a ogr dataset (most supported: shp, geojson, gml, etc)
    else if (strcmp(argv[argNum], "--ogr") == 0) {
#if GDAL_VERSION_MAJOR < 2
      OGRRegisterAll();
      OGRDataSource *dataSource = OGRSFDriverRegistrar::Open(argv[argNum+1], false);
#else
      GDALAllRegister();
      GDALDataset *dataSource = (GDALDataset*) GDALOpenEx(argv[argNum+1], GDAL_OF_READONLY, NULL, NULL, NULL);
#endif
      ++argNum;
      if (dataSource == NULL) {
        std::cerr << "Error: Could not open file." << std::endl;
        return false;
      }
      OGRLayer *dataLayer = dataSource->GetLayer(0); //-- get first layer
      dataLayer->ResetReading();
      //-- Reads all features in this layer
      OGRFeature *feature;
      //-- get first feature
      if (dataLayer->GetFeatureCount(true) > 1)
        std::cout << "Reading only the first feature in the shapefile." << std::endl;
      feature = dataLayer->GetNextFeature();
      if (feature->GetGeometryRef()->getGeometryType() == wkbPolygon) {
        geometry = static_cast<OGRPolygon *>(feature->GetGeometryRef());feature->GetGeometryRef();
      }
      else {
        std::cout << "First feature ain't a POLYGON." << std::endl;
        return(0);
      }
      
    }
    
    else {
      usage();
      return(0);
    }
  }
  
  if (!timeResults) startTime = 0;
  
  PolygonRepair prepair;
  
  //-- compute robustness
  if (computeRobustness == true)
    std::cout << "Robustness of input polygon: " << sqrt(prepair.computeRobustness(geometry)) <<std::endl;
  
  
  OGRMultiPolygon *outPolygons;
  if (pointSet) {
    outPolygons = prepair.repairPointSet(geometry, startTime);
  } else {
    outPolygons = prepair.repairOddEven(geometry, startTime);
  }
  
  if (minArea > 0) {
    prepair.removeSmallPolygons(outPolygons, minArea);
  }
  
  //-- output well known text
  if (shpOut) {
    prepair.saveToShp(outPolygons, "out.shp");
  }
  else {
    char *outputWKT;
    outPolygons->exportToWkt(&outputWKT);
    std::cout << outputWKT << std::endl;
  }
  
  
  //-- compute robustness
  if (computeRobustness == true)
    std::cout << "Robustness of output polygon: " << sqrt(prepair.computeRobustness()) <<std::endl;
  
  //-- time results
  if (timeResults) {
    time_t totalTime = time(NULL)-startTime;
    std::cout << "Done! Process finished in " << totalTime/60 << " minutes " << totalTime%60 << " seconds." << std::endl;
  }
  
  return 0;
}

void usage() {
  std::cout << "=== prepair Help ===\n" << std::endl;
  std::cout << "Usage:   prepair --wkt 'POLYGON((0 0, 0 10, 10 0, 10 10, 0 0))'" << std::endl;
  std::cout << "OR" << std::endl;
  std::cout << "Usage:   prepair -f infile.txt (infile.txt must contain one WKT on the 1st line)" << std::endl;
  std::cout << "OR" << std::endl;
  std::cout << "Usage:   prepair --ogr infile.shp (first polygon of infile.shp is processed)" << std::endl;
  std::cout << "================================================================================" << std::endl;
  std::cout << "Additional options:" << std::endl;
  std::cout << "--robustness   Compute the robustness of the input/output" << std::endl;
  std::cout << "--setdiff      Uses the point set topology paradigm (default: odd-even paradigm)" << std::endl;
  std::cout << "--minarea AREA Only output polygons larger than AREA" << std::endl;
  std::cout << "--isr GRIDSIZE Snap round the input before repairing" << std::endl;
  std::cout << "--time         Benchmark the different stages of the process" << std::endl;
  std::cout << "--shpOut       Output a shapefile (out.shp) instead of a WKT" << std::endl;
  
}


