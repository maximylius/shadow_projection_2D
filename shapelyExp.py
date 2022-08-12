# -*- coding: utf-8 -*-
"""
Script to create Sunlight / Shadow Maps


"""

import geopandas as gpd
import matplotlib.pyplot as plt
import shapely
from shapely.geometry import Polygon, LineString
from pysolar import solar
import datetime
import math
import pytz









from shapely.geometry import  MultiLineString, mapping, shape
coords = [((0, 0), (1, 1)), ((-1, 0), (1, 0))]
coords = [((0, 0), (1, 1)), ((-1, 0), (1, 0))]
coords = [(0, 0), (1, 1), (-1, 0), (1, 0)]
# geom.bounds
lines = MultiLineString(coords)
lines 
print(lines)
for line in lines:
    line
    print(line)

# convert to GeoJSON format:
print(mapping(lines))
{'type': 'MultiLineString', 'coordinates': (((0.0, 0.0), (1.0, 1.0)), ((-1.0, 0.0), (1.0, 0.0)))}
# convert from GeoJSON to shapely
print(shape(mapping(lines)))
# MULTILINESTRING ((0 0, 1 1), (-1 0, 1 0))



c0 = [(0,0), (4,0), (4,1), (1,1), (1,3), (0,3)]
s0 = Polygon(c0)    
s0
s1 = Polygon(map(lambda x: (x[0]+3, x[1]+0) , c0))
shapely.ops.cascaded_union([s0,s1])    
# shapely.ops.cascaded_union([s0,s1]).convex_hull
l0 = LineString([(0,0),(2,2)])
l1 = LineString([(1,1),(3,3)])
l1
list(shapely.affinity.scale(l1, 3, 3).coords)
# l2 = LineString([(2,2),(4,4)])
# l3 = LineString([(3,3),(5,5)])
# l4 = LineString([(2,0),(0,2)])
# l0.covers(l1)
# l0.intersects(l1)
# l0.overlaps(l1)
# l0.crosses(l1)
# l0.contains(l1)
# l0.disjoint(l1)
# l0.disjoint(l2)
# l0.disjoint(l3)
# d01 = l0.difference(l1)    
# l0.difference(l2) == l0 # difference equal original -> no overlap
# list(d01.coords)
# sd01 = l0.symmetric_difference(l1)
# l0.touches(l2) # true if they only touch in one point: intersection but not crosses / overlap
# shapely.ops.cascaded_union([geom, geomShifted]).convex_hull
# 2do: define REAL intersections correctly
from shapely.affinity import translate 
sq = Polygon([(0,0), (0,1), (1,1), (1,0)])
sqU = shapely.ops.cascaded_union([
    translate(sq,0,0),
    translate(sq,1,0),
    translate(sq,2,0),
    translate(sq,2,1),
    translate(sq,0,1)
    ])

sqO = shapely.ops.cascaded_union([
    translate(sq,0,0),
    translate(sq,1,0),
    translate(sq,2,0),
    translate(sq,2,1),
    translate(sq,0,1),
    translate(sq,0,2),
    translate(sq,1,2),
    translate(sq,2,2)
    ])

sqj = shapely.ops.cascaded_union([
    translate(sq,0,0),
    translate(sq,1,0),
    translate(sq,2,0),
    translate(sq,2,1),
    translate(sq,2,2)
    ])
xoff=0.1
yoff=0.2
nSteps = 11
sqjP = shapely.ops.cascaded_union(
    list(map(lambda x: translate(sqj, x*xoff, x*yoff), range(nSteps)))
)
sqjP 
sqjP.simplify((xoff**2+yoff**2)**0.5)

sqjT = Polygon([
    (0,0),
    (0,1),
    (0+nSteps*xoff, 1+nSteps*yoff),
    (0+nSteps*xoff, 1+nSteps*yoff),
    
    nSteps
    ])

sqU.simplify(0.5)
sqO.simplify(20)

# Polygon([[0,0],[2,0], [2,2], [1,2], [1,1],[1,2],  [0,2]]) # invalid
polygon = Polygon([[1,0],[2,5], [4,4], [4,0]])
path = LineString([(3,3), (5,5)])
path.intersects(polygon)
pathInt = list(path.intersection(polygon).coords)
comb = list(path.coords)+pathInt
comb.sort()
comb


import geopandas as gpd
from shapely.geometry import shape, MultiLineString
import pandas as pd

import random
random.uniform(-10, 10)

n = 40
randomLines = list(map(
    lambda x: LineString([
        (random.uniform(-10, 10),random.uniform(-10, 10)),
        (random.uniform(-10, 10),random.uniform(-10, 10))
        ]), 
    range(n)
    ))
gpd.GeoSeries(randomLines ).plot()


lines = shapely.geometry.MultiLineString(randomLines)
## Convert Lines to Polygons by applying a tiny buffer
buf = 0.0000000000001
lines = lines.buffer(buf)
## Get outer boundary of the lines as a polygon
boundary = lines.convex_hull
## Get a difference to generate a multipolygon
multipolygons = list(boundary.difference(lines))
wb = list(map(
    lambda x: x.buffer(buf*2),
    list(multipolygons)
    ))
plotGeometry(list(wb), None, '', 1)

# geom = polygon

# path = LineString([(1,0), (4,3)])
# path2 = LineString([(2,1), (3,2)])
# path.intersects(path2)

# rotation not necessary to find left, right and bottom
#       calculate signed difference of points to line
#       distance parallel to azimuth angle for left and right most
#       distance to orthogonal line to azimuth angle for bottom point