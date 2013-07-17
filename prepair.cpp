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

#include "PolygonRepair.h"

//-- minimum size of a polygon in the output (smaller ones are not returned)
//-- can be changed with --min flag
double MIN_AREA = 0;
double ISR_TOLERANCE = 0;
bool wktout = false;
bool shpout = false;

void usage();
bool savetoshp(OGRMultiPolygon* multiPolygon);

int main (int argc, const char * argv[]) {
    
    if (argc < 2 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        usage();
        return(0);
    }
    
    OGRGeometry *geometry;
    
    for (int argNum = 1; argNum < argc; ++argNum) {
        if (strcmp(argv[argNum], "--minarea") == 0) {
            MIN_AREA = atof(argv[argNum+1]);
            ++argNum;
        }
        //-- ISR snapping tolerance
        else if (strcmp(argv[argNum], "--isr") == 0) {
            ISR_TOLERANCE = atof(argv[argNum+1]);
            ++argNum;
            // TODO : scale dataset if tolerance < 1.0
        }
        //-- reading from WKT passed directly
        else if (strcmp(argv[argNum], "--wkt") == 0) {
            unsigned int bufferSize = 100000000;
            char *inputWKT = (char *)malloc(bufferSize*sizeof(char *));
            strcpy(inputWKT, argv[argNum+1]);
            ++argNum;
            OGRGeometryFactory::createFromWkt(&inputWKT, NULL, &geometry);
            if (geometry == NULL) {
                std::cout << "Error: WKT is not valid" << std::endl;
                return 1;
            } wktout = true;
        }
        //-- reading from WKT stored in first line of a text file
        else if (strcmp(argv[argNum], "-f") == 0) {
            unsigned int bufferSize = 100000000;
            char *inputWKT = (char *)malloc(bufferSize*sizeof(char *));
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
            } wktout = true;
        }
        //-- reading from a shapefile
        else if (strcmp(argv[argNum], "--shp") == 0) {
            OGRRegisterAll();
            OGRDataSource *dataSource = OGRSFDriverRegistrar::Open(argv[argNum+1], false);
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
                shpout = true;
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
    
    PolygonRepair prepair;
    OGRMultiPolygon *outPolygons = prepair.repairOddEven(geometry);
    
    //-- snap rounding of the input
    /*if (ISR_TOLERANCE != 0) {
        Polyline_list_2 *snappedgeom = isr(geometry);
        std::cout << "done." << std::endl;
        outPolygons = repair(snappedgeom);
    }
    else {
        outPolygons = repair(geometry);
    }*/
    
    
    char *outputWKT;
    /*if (MIN_AREA > 0) {
        std::cout << "Removing polygons smaller than " << MIN_AREA << " unit^2." << std::endl;
        for (std::list<OGRPolygon*>::iterator it = outPolygons->begin(); it != outPolygons->end(); ++it) {
            if ((*it)->get_Area() > MIN_AREA) {
                multiPolygon->addGeometryDirectly(*it);
            }
        }
    }
    else {
        for (std::list<OGRPolygon*>::iterator it = outPolygons->begin(); it != outPolygons->end(); ++it)
            multiPolygon->addGeometryDirectly(*it);
    }*/
    
    if (wktout) {
        outPolygons->exportToWkt(&outputWKT);
        std::cout << std::endl << "Repaired polygon:" << std::endl << outputWKT << std::endl;
    }
    
    //-- save to a shapefile
    else if (shpout) {
        savetoshp(outPolygons);
    } return 0;
}

void usage() {
    std::cout << "=== prepair Help ===\n" << std::endl;
    std::cout << "Usage:   prepair --wkt 'POLYGON(...)'" << std::endl;
    std::cout << "OR" << std::endl;
    std::cout << "Usage:   prepair -f infile.txt (infile.txt must contain one WKT on the 1st line)" << std::endl;
    std::cout << "OR" << std::endl;
    std::cout << "Usage:   prepair --shp infile.shp (first polygon of infile.shp is processed)" << std::endl;
    
}


bool savetoshp(OGRMultiPolygon* multiPolygon) {
    const char *driverName = "ESRI Shapefile";
    OGRRegisterAll();
	OGRSFDriver *driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(driverName);
	if (driver == NULL) {
		std::cout << "\tError: OGR Shapefile driver not found." << std::endl;
		return false;
	}
	OGRDataSource *dataSource = driver->Open("out.shp", false);
	if (dataSource != NULL) {
		std::cout << "\tOverwriting file..." << std::endl;
		if (driver->DeleteDataSource(dataSource->GetName())!= OGRERR_NONE) {
			std::cout << "\tError: Couldn't erase file with same name." << std::endl;
			return false;
		} OGRDataSource::DestroyDataSource(dataSource);
	}
	std::cout << "\tWriting file... " << std::endl;
	dataSource = driver->CreateDataSource("out.shp", NULL);
	if (dataSource == NULL) {
		std::cout << "\tError: Could not create file." << std::endl;
		return false;
	}
	OGRLayer *layer = dataSource->CreateLayer("polygons", NULL, wkbPolygon, NULL);
	if (layer == NULL) {
		std::cout << "\tError: Could not create layer." << std::endl;
		return false;
	}
    OGRFeature *feature = OGRFeature::CreateFeature(layer->GetLayerDefn());
    // feature->SetField("Name", szName);
    feature->SetGeometry(multiPolygon);
    if (layer->CreateFeature(feature) != OGRERR_NONE) {
        std::cout << "\tError: Could not create feature." << std::endl;
    }
    OGRFeature::DestroyFeature(feature);
    OGRDataSource::DestroyDataSource(dataSource);
    return true;
}


/*Polyline_list_2* isr(OGRGeometry* geometry) {
 std::cout << "ISR snapping with tolerance: " << ISR_TOLERANCE << std::endl;
 Segment_list_2 seg_list;
 OGRPolygon *polygon = (OGRPolygon *)geometry;
 polygon->closeRings();
 for (int currentPoint = 0; currentPoint < (polygon->getExteriorRing()->getNumPoints() - 1); ++currentPoint) {
 Segment_2 s = Segment_2( ISRPoint(polygon->getExteriorRing()->getX(currentPoint), polygon->getExteriorRing()->getY(currentPoint)),
 ISRPoint(polygon->getExteriorRing()->getX(currentPoint+1), polygon->getExteriorRing()->getY(currentPoint+1)));
 if (s.is_degenerate() == false)
 seg_list.push_back(s);
 // else
 // std::cout << "degenerate segment" << std::endl;
 }
 // ISRPoint(polygon->getExteriorRing()->getX((currentPoint+1)%polygon->getExteriorRing()->getNumPoints()),
 //       polygon->getExteriorRing()->getY((currentPoint+1)%polygon->getExteriorRing()->getNumPoints()))
 // TODO: @Ken: I think you're creating twice the same edge here... no?
 
 
 for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing) {
 for (int currentPoint = 0; currentPoint < (polygon->getInteriorRing(currentRing)->getNumPoints() - 1); ++currentPoint) {
 Segment_2 s = Segment_2(ISRPoint(polygon->getInteriorRing(currentRing)->getX(currentPoint),   polygon->getInteriorRing(currentRing)->getY(currentPoint)),
 ISRPoint(polygon->getInteriorRing(currentRing)->getX(currentPoint+1), polygon->getInteriorRing(currentRing)->getY(currentPoint+1)));
 if (s.is_degenerate() == false)
 seg_list.push_back(s);
 // else
 // std::cout << "Degenerate segment skipped." << std::endl;
 }
 }
 
 // std::cout << "input segments: " << seg_list.size() << std::endl;
 // Segment_list_2::const_iterator its;
 // for (its = seg_list.begin(); its != seg_list.end(); ++its) {
 //   std::cout << "Segment " << (*its) << std::endl;
 //   std::cout << "   " << sqrt(CGAL::to_double(its->squared_length())) << std::endl;
 // }
 
 Polyline_list_2 *output_list = new Polyline_list_2();
 // Polyline_list_2 output_list;
 CGAL::snap_rounding_2<Traits,Segment_list_2::const_iterator, Polyline_list_2>
 (seg_list.begin(), seg_list.end(), *output_list, ISR_TOLERANCE, true, false, 2);
 
 int counter = 0;
 Polyline_list_2::const_iterator iter1;
 for (iter1 = output_list->begin(); iter1 != output_list->end(); ++iter1) {
 std::cout << "Polyline number " << ++counter << ":\n";
 Polyline_2::const_iterator iter2;
 // std::cout << iter1->size() << std::endl;
 for (iter2 = iter1->begin(); iter2 != iter1->end(); ++iter2)
 std::cout << CGAL::to_double(iter2->x()) << " " << CGAL::to_double(iter2->y()) << std::endl;
 }
 
 return output_list;
 // return NULL;
 }*/


/*std::list<OGRPolygon*>* repair(Polyline_list_2* seg_list) {
 // Triangulation
 Triangulation triangulation;
 
 // int counter = 0;
 Polyline_list_2::const_iterator iter1;
 for (iter1 = seg_list->begin(); iter1 != seg_list->end(); iter1++) {
 // std::cout << "Polyline number " << ++counter << ":\n";
 Polyline_2::const_iterator iter2;
 ISRPoint last = iter1->back();
 for (iter2 = iter1->begin(); iter2 != iter1->end(); ) {
 if (*iter2 == last)
 break;
 Point p1 = Point(CGAL::to_double(iter2->x()), CGAL::to_double(iter2->y()));
 ++iter2;
 Point p2 = Point(CGAL::to_double(iter2->x()), CGAL::to_double(iter2->y()));
 triangulation.insert_constraint(p1, p2);
 }
 }
 
 std::list<OGRPolygon*>* outPolygons = repair_tag_triangulation(triangulation);
 return outPolygons;
 }*/