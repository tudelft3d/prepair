## What is prepair?

prepair — pronounce 'pee-repair' as in 'polygon repair' — permits us to easily repair "broken" GIS polygons according to the international standard [ISO19107](http://www.iso.org/iso/catalogue_detail.htm?csnumber=26012) (Geographic information — Spatial schema). Given one input polygon, it *automatically* repairs it and returns back a valid polygon (actually a MultiPolygon since the input can represent more than one polygon — think of a 'bowtie' for instance). 

Automated repair methods can be considered as interpreting ambiguous or ill-defined polygons and giving a coherent and clearly defined output. Examples of errors are: polygon has a dangling edge; polygon is not closed; polygon self-intersects; an inner ring of the polygon is located outside the outer ring; etc.

prepair performs more or less the same as the PostGIS 2.0's function [ST_MakeValid()](http://postgis.org/documentation/manual-svn/ST_MakeValid.html), but is faster, scales better to massive polygons, and predicting its behaviour is simple (so one can guess how her polygons will be repaired).
We have implemented two repair paradigms: 

  1. an extension of the [odd-even algorithm](https://en.wikipedia.org/wiki/Even-odd_rule) to handle GIS polygons containing inner rings and degeneracies; 
  2. setdiff: one where we follow a *point set difference* rule for the rings (outer - inner).

prepair is based on a constrained triangulation ([CGAL](http://www.cgal.org) is used) and [OGR](http://www.gdal.org/ogr/) is used to read/write WKT.

It is available under a dual license scheme: [GPLv3](http://www.gnu.org/copyleft/gpl.html) and commercial. If you are interested in a commercial license, please contact [Ken Arroyo Ohori](mailto:g.a.k.arroyoohori@tudelft.nl).

Note that prepair is only concerned with single polygons, and if you're interested in validating how different polygons interact with each other (to be precise: to check if they form a planar partition) have a look at our other project [pprepair](https://github.com/tudelft3d/pprepair).


## Details

Details of how we automatically repair broken polygons, and what results you can expect, are available in this scientific article:

> Ledoux, H., Arroyo Ohori, K., and Meijers, M. (2014). A triangulation-based approach to automatically repair GIS polygons. *Computers & Geosciences* 66:121–131. [ [DOI] ](http://dx.doi.org/10.1016/j.cageo.2014.01.009) [ [PDF] ](http://3dgeoinfo.bk.tudelft.nl/hledoux/pdfs/14_cgeo_prepair.pdf)

If you use prepair for a scientific project, please cite this article.


## How to get it?

prepair is provided as source code or as 64-bit binaries for [Windows](https://github.com/tudelft-gist/prepair/releases/download/v0.7/prepair_win64.zip) and [Mac](https://github.com/tudelft-gist/prepair/releases/download/v0.7/prepair_mac.zip). The Mac binary requires Kyngchaos' [GDAL 1.11 Complete Framework](http://www.kyngchaos.com/software/frameworks#gdal_complete).

prepair is also very easy to compile on Mac and Linux using the included CMake file. It should also work on other Unix-like systems. To compile prepair, you need to have a recent version of the following three (free) libraries:

1. [CGAL](http://www.cgal.org)
2. [OGR](http://www.gdal.org/ogr/)
3. [CMake](http://www.cmake.org) 

Under Mac, if you use Kyngchaos' GDAL Complete Framework, which is used by QGIS, you already have OGR installed. If you need them, a good way to install CGAL and OGR is to use [Homebrew](http://brew.sh):

    $ brew install gdal
    $ brew install cgal 

Once all the dependencies are met, just generate the makefile for your system and compile:

    $ cmake .
    $ make


## How to run it?

You can run prepair from the command-line or through our [QGIS plug-in](https://github.com/tudelft-gist/prepair-qgis), which you can get from the official QGIS repository. 

A [WKT](http://en.wikipedia.org/wiki/Well-known_text) or an OGR dataset (shapefile, geojson or GML for instance) is read as input, and a WKT or a shapefile (a MultiPolygon) is given as output:

    $ ./prepair --wkt 'POLYGON((0 0, 0 10, 10 0, 10 10, 0 0))'  
    MULTIPOLYGON (((0 10,0 0,5 5,0 10)),((5 5,10 0,10 10,5 5)))  
    
    $ ./prepair --ogr myfile.shp
    MULTIPOLYGON (((0 10,0 0,5 5,0 10)),((5 5,10 0,10 10,5 5)))

    $ ./prepair --shpOut --ogr data/CLC2006_180927.geojson 
    Creating out.shp
    

[Snap rounding](http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Snap_rounding_2/Chapter_main.html) of the input segments can be performed with the --isr option:

    $ ./prepair --isr 2 --wkt "POLYGON((0 0, 10 0, 15 5, 10 0, 10 10, 0 10, 0 0))"
    MULTIPOLYGON (((11 1,11 11,1 11,1 1,11 1)))
    

It's also possible to remove small (sliver) polygons in the output by giving the smallest area allowed with the --minarea option:

    $ ./prepair --wkt 'POLYGON((0 0, 10 0, 10 11, 11 10, 0 10))' 
    MULTIPOLYGON (((10 0,10 10,0 10,0 0,10 0)),((11 10,10 11,10 10,11 10)))

    $ ./prepair --wkt 'POLYGON((0 0, 10 0, 10 11, 11 10, 0 10))' --minarea 1
    Removing polygons smaller than 1 unit^2.
    MULTIPOLYGON (((10 0,10 10,0 10,0 0,10 0)))


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

