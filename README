## What is prepair?

prepair permits you to easily repair "broken" GIS polygons, and that according to the international standards ISO 19107. In brief, given a polygon stored in [WKT](http://en.wikipedia.org/wiki/Well-known_text), it *automatically* repairs it and gives you back a valid WKT. Automated repair methods can be considered as interpreting ambiguous or ill-defined polygons and giving a coherent and clearly defined output.

It performs more or less the same as the new PostGIS 2.0's function [ST_MakeValid()](http://postgis.org/documentation/manual-svn/ST_MakeValid.html), but is several order of magnitude faster, scales better to massive polygons, and predicting its behaviour is simple (so one can guess how her polygons will be repaired).

prepair is based on a constrained triangulation ([CGAL](http://www.cgal.org) is used) and [OGC](http://www.gdal.org/ogr/) is used to read/write WKT.

## Details
Details of how we automatically repair broken polygons, and what results you can expect, are available our [Agile 2012 paper](http://www.gdmc.nl/ledoux/pdfs/_12agile.pdf).

## How to compile?

You need to install the following two free libraries:

1. [CGAL](http://www.cgal.org)
2. [OGR](http://www.gdal.org/ogr/)

And then use the makefile provided for Mac and Linux. For Windows, you're on your own right now, but we plan to provide binaries in the near future.

## It's a command-line interface only

WKT are read as input, and a WKT (a MultiPolygon) is given as output:

    $ ./prepair 'POLYGON((0 0, 0 10, 10 0, 10 10, 0 0))'  
    Processing: POLYGON((0 0, 0 10, 10 0, 10 10, 0 0))  
    Repaired polygon:  
    MULTIPOLYGON (((0 10,0 0,5 5,0 10)),((5 5,10 0,10 10,5 5)))  
    
## Examples of input you can try ##

A 'bowtie' polygon: POLYGON((0 0, 0 10, 10 0, 10 10, 0 0))

Square with wrong orientation: POLYGON((0 0, 0 10, 10 10, 10 0, 0 0))



