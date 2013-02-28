## What is prepair?

prepair -- pronounce 'pee-repair' as in 'polygon repair' -- permits us to easily repair "broken" GIS polygons, and that according to the international standards ISO 19107 (Geographic information -- Spatial schema). In brief, given a polygon stored in [WKT](http://en.wikipedia.org/wiki/Well-known_text), it *automatically* repairs it and returns back a valid WKT. Note that this project is only concerned with single polygons, if you're interested in validating how different polygons interact with each others (to be precise: to check if they form a planar partition) have a look at our other project [pprepair](https://github.com/tudelft-gist/pprepair).

Automated repair methods can be considered as interpreting ambiguous or ill-defined polygons and giving a coherent and clearly defined output. Examples of errors are: polygon has a dangling edge; polygon is not closed; polygon self-intersects; an inner ring of the polygon is located outside the outer ring; etc.

It performs more or less the same as the new PostGIS 2.0's function [ST_MakeValid()](http://postgis.org/documentation/manual-svn/ST_MakeValid.html), but is several order of magnitude faster, scales better to massive polygons, and predicting its behaviour is simple (so one can guess how her polygons will be repaired).

prepair is based on a constrained triangulation ([CGAL](http://www.cgal.org) is used) and [OGR](http://www.gdal.org/ogr/) is used to read/write WKT.

## Details
Details of how we automatically repair broken polygons, and what results you can expect, are available in our [Agile 2012 paper](http://www.gdmc.nl/ledoux/pdfs/_12agile.pdf).

## How to compile?

You first need to install the following two (free) libraries:

1. [CGAL](http://www.cgal.org)
2. [OGR](http://www.gdal.org/ogr/)

And then use the makefile provided for Mac and Linux. For Windows, you're on your own right now, but we plan to provide binaries in the near future.

## It's a command-line program only

WKT is read as input, and a WKT (a MultiPolygon) is given as output:

    $ ./prepair 'POLYGON((0 0, 0 10, 10 0, 10 10, 0 0))'  
    Processing: POLYGON((0 0, 0 10, 10 0, 10 10, 0 0))  
    Repaired polygon:  
    MULTIPOLYGON (((0 10,0 0,5 5,0 10)),((5 5,10 0,10 10,5 5)))  
    
## Examples of invalid input you can try ##

A 'bowtie' polygon: 
    
    POLYGON((0 0, 0 10, 10 0, 10 10, 0 0))

Square with wrong orientation: 
    
    POLYGON((0 0, 0 10, 10 10, 10 0, 0 0))

Inner ring touching the outer ring:

    POLYGON((0 0, 10 0, 10 10, 0 10, 0 0),(5 2,5 7,10 7, 10 2, 5 2))

Dangling edge:

    POLYGON((0 0, 10 0, 15 5, 10 0, 10 10, 0 10, 0 0))

Outer ring not closed:

    POLYGON((0 0, 10 0, 10 10, 0 10))

Two adjacent inner rings:

    POLYGON((0 0, 10 0, 10 10, 0 10, 0 0), (1 1, 1 8, 3 8, 3 1, 1 1), (3 1, 3 8, 5 8, 5 1, 3 1))

Polygon with an inner ring inside another inner ring:

    POLYGON((0 0, 10 0, 10 10, 0 10, 0 0), (2 8, 5 8, 5 2, 2 2, 2 8), (3 3, 4 3, 3 4, 3 3))

