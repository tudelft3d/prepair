## What is prepair?

prepair — pronounce 'pee-repair' as in 'polygon repair' — permits us to easily repair "broken" GIS polygons according to the international standard [ISO19107](http://www.iso.org/iso/catalogue_detail.htm?csnumber=26012) (Geographic information — Spatial schema). Given one input polygon, it *automatically* repairs it and returns back a valid polygon (actually a MultiPolygon since the input can represent more than one polygon — think of a 'bowtie' for instance). 

Automated repair methods can be considered as interpreting ambiguous or ill-defined polygons and giving a coherent and clearly defined output. Examples of errors are: polygon has a dangling edge; polygon is not closed; polygon self-intersects; an inner ring of the polygon is located outside the outer ring; etc.

The algorithm behind prepair was integrated into [CGAL](https://www.cgal.org) as the [2D Polygon Repair](https://doc.cgal.org/latest/Polygon_repair/index.html) package (CGAL 6.0+). This software is now a thin wrapper around that package, adding GDAL/OGR I/O support and 3D polygon repair via plane fitting. For direct use in C++ projects, consider using the CGAL package directly.

prepair uses [GDAL](http://www.gdal.org/) to read/write geometries and [Boost](http://www.boost.org/) for command-line parsing.

It is available under the [GPLv3](http://www.gnu.org/copyleft/gpl.html) licence, which allows you to use, copy and modify the software freely. However, if you incorporate prepair in your software, you must distribute the source code of your software, as well as any modifications made to pprepair, under the GPLv3 as well.

Note that prepair is only concerned with single polygons, and if you're interested in validating how different polygons interact with each other (to be precise: to check if they form a planar partition) have a look at our other project [pprepair](https://github.com/tudelft3d/pprepair).

If you want to use the polygon repair algorithm directly in your C++ code, see the [CGAL 2D Polygon Repair](https://doc.cgal.org/latest/Polygon_repair/index.html) package instead.


## Details

Details of how we automatically repair broken polygons, and what results you can expect, are available in this scientific article:

> Ledoux, H., Arroyo Ohori, K., and Meijers, M. (2014). A triangulation-based approach to automatically repair GIS polygons. *Computers & Geosciences* 66:121–131. [ [DOI] ](http://dx.doi.org/10.1016/j.cageo.2014.01.009) [ [PDF] ](http://3dgeoinfo.bk.tudelft.nl/hledoux/pdfs/14_cgeo_prepair.pdf)

If you use prepair for a scientific project, please cite this article.

## How to get it?

prepair is very easy to compile on Mac and Linux using the included CMake file. It should also work on other Unix-like systems and is possible to compile under Windows. To compile prepair, you need the following libraries:

1. [CGAL](http://www.cgal.org) 6.0 or later (provides the polygon repair algorithm)
2. [GDAL](http://www.gdal.org/) 3.0 or later (for I/O)
3. [Boost](http://www.boost.org/) 1.50 or later (for command-line parsing)
4. [CMake](http://www.cmake.org) 

On Mac, you can install it using [Homebrew](http://brew.sh):

    $ brew install tudelft3d/software/prepair

Once all the dependencies are met, just generate the makefile for your system and compile:

    $ cmake -DCMAKE_BUILD_TYPE=Release .
    $ make

## How to run it?

You can run prepair from the command-line or through our [QGIS plug-in](https://github.com/tudelft-gist/prepair-qgis), which you can get from the official QGIS repository. 

A [WKT](http://en.wikipedia.org/wiki/Well-known_text) or a path to a dataset (geopackage, geojson or shapefile for instance) is read as input, and a WKT or a path to a dataset is given as output:

    $ ./prepair -w 'POLYGON((0 0, 0 10, 10 0, 10 10, 0 0))'  
    MULTIPOLYGON (((0 0,5 5,0 10,0 0)),((5 5,10 0,10 10,5 5)))
    
    $ ./prepair -w 'POLYGON((0 0, 0 10, 10 0, 10 10, 0 0))' -o myfile.gpkg
    Writing GPKG file myfile.gpkg...
    
    $ ./prepair -w 'POLYGON((0 0 0, 0 10 5, 10 0 10, 10 10 15, 0 0 0))'
    MULTIPOLYGON (((0 0,5 5,0 10,0 0)),((5 5,10 0,10 10,5 5)))

## Examples of invalid input you can try

The folder 'data' contains examples of relatively big invalid polygons. These are from the [Corine Land Cover 2006 dataset](http://sia.eionet.europa.eu/CLC2006).

A 'bowtie' polygon: 
    
    POLYGON((0 0, 0 10, 10 0, 10 10, 0 0))

Square with wrong orientation: 
    
    POLYGON((0 0, 0 10, 10 10, 10 0, 0 0))

Inner ring with one edge sharing part of an edge of the outer ring:

    POLYGON((0 0, 10 0, 10 10, 0 10, 0 0),(5 2,5 7,10 7, 10 2, 5 2))

Dangling edge:

    POLYGON((0 0, 10 0, 15 5, 10 0, 10 10, 0 10, 0 0))

Outer ring not closed:

    POLYGON((0 0, 10 0, 10 10, 0 10))

Two adjacent inner rings:

    POLYGON((0 0, 10 0, 10 10, 0 10, 0 0), (1 1, 1 8, 3 8, 3 1, 1 1), (3 1, 3 8, 5 8, 5 1, 3 1))

Polygon with an inner ring inside another inner ring:

    POLYGON((0 0, 10 0, 10 10, 0 10, 0 0), (2 8, 5 8, 5 2, 2 2, 2 8), (3 3, 4 3, 3 4, 3 3))

