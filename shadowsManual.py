# -*- coding: utf-8 -*-
"""
Script to create Sunlight / Shadow Maps


"""


# cd /downloads
# pip install GDAL-3.3.1-cp38-cp38-win_amd64.whl
# pip install pyproj-3.1.0-cp38-cp38-win_amd64.whl
# pip install Fiona-1.8.20-cp38-cp38-win_amd64.whl
# pip install Shapely-1.7.1-cp38-cp38-win_amd64.whl
# pip install geopandas-0.9.0-py3-none-any.whl


# %% set up
# set path
import os
#path = 'C:\\Users\\BSE\\Documents\\MaxVonMylius\\GeographyGrowth\\TermPaper\\py'
path = 'C:\\Users\\mmyli\\Documents\\HumboldtUni\\GeographyGrowth\\Paper\\shadows'
os.chdir(path)

import geopandas as gpd
import matplotlib.pyplot as plt
import shapely
from shapely.geometry import Polygon, LineString



# %% load shapefile

# Importing an ESRI Shapefile
# veg = gpd.read_file(r'./../Sample/sampleVeg.shp')
# geb = gpd.read_file(r'./../Sample/sampleGeb.shp')
veg = gpd.read_file(r'./Sample/sampleVeg.shp')
geb = gpd.read_file(r'./Sample/sampleGeb.shp')


# check projection
print(veg.crs) # EPSG: 25833
print(geb.crs) # EPSG: 25833

# reprojecting to DHDN / Soldner Berlin EPSG:3068 with transformation: 1777
crs_code = 3068
veg = veg.to_crs(epsg = crs_code)
geb = geb.to_crs(epsg = crs_code)



# plot using GeoPandas and Matplotlib
fig, axs = plt.subplots(nrows = 2, ncols = 2, figsize = (12,8))
veg.plot(ax=axs[0,0])
veg.plot(ax=axs[0,1], cmap = 'Greens', edgecolor = 'black', column = 'Mean_nDOM')
geb.plot(ax=axs[1,0], cmap = 'Reds', edgecolor = 'black', column = 'Mean_nDOM')
geb.plot(ax=axs[1,1], cmap = 'Reds', edgecolor = 'black', column = 'Max_nDOM')
veg.plot(ax=axs[0,0], color = "#11442288", edgecolor = "pink")


geb
geb.iloc[20, list(geb.columns).index('geometry')]

# %% set up sunlight angles
from pysolar import solar
import datetime
import math
import pytz

# year, month, day, hour, minute, second, microsecond
date = datetime.datetime(2007, 2, 16, 12, 13, 1, 130320, tzinfo=pytz.timezone('Europe/Berlin'))


latitude_southMostBerlin = 52.388447
longitude_westMostBerlin = 13.131997
latitude_northMostBerlin = 52.675258
longitude_northMostBerlin = 13.479466
print(solar.get_azimuth(latitude_southMostBerlin, longitude_westMostBerlin, date))
print(solar.get_azimuth(latitude_northMostBerlin, longitude_northMostBerlin, date))

print(solar.get_altitude(latitude_southMostBerlin, longitude_westMostBerlin, date))
print(solar.get_altitude(latitude_northMostBerlin, longitude_northMostBerlin, date))


latitude = 42.206
longitude = -71.382



print(date)
# calculate angle between sun and earth surface at location and date
altitude = solar.get_altitude(latitude, longitude, date) 
# calculate direction of sunlight at location and date
azimuth = solar.get_azimuth(latitude, longitude, date)
print(altitude)
print(azimuth)

# if azimuth between 0° and 180° positive (=going east), else negative
x = math.sin((azimuth) * math.pi / 180)
# if azimuth between 270° and 90° positive (=going north) else negative
y = math.sin((azimuth - 270) * math.pi / 180)

offsetMultiplicator = math.tan(altitude * math.pi / 180)
cardinalDirection = [x,y] #[north=1-south=-1, west=1-east=-1,]


# %% create projections

vegShadows = veg[
    (veg['Mean_nDOM'] > 0) & 
    (veg['Class_name'] != 'Veg < 2,5 m (einschl. Rasen, Grünland)')
]

vegShadows['height'] = vegShadows['Mean_nDOM']

#reproject to obtain proper WGS84 latitude-longitude projection EPSG:4326, as proper
vegShadows = vegShadows.to_crs(epsg = 4326)
vegShadows['centroid'] = vegShadows.centroid
vegShadows['long4326'] = vegShadows['centroid'].x
vegShadows['lat4326'] = vegShadows['centroid'].y
# project back to berlin projection
vegShadows = vegShadows.to_crs(epsg = crs_code)

# vegShadow = vegShadow[vegShadow.index < vegShadow.index[2]]

date = datetime.datetime(2007, 4, 16, 16, 13, 1, 130320, tzinfo=pytz.timezone('Europe/Berlin'))

# helper functions outside scope
def getEdgesFromCoords (coords):    
    edges = []
    # 2do: Problem! there might be groups of edges (hole in the building) 
    # thus not all points shall be connected
    # only points within groups of edges
    # len(list(geom.interiors)) == 0 # no holes
    for c0 in range(len(coords)):
        # define index + 1
        c1 = c0 + 1 if c0 + 1 < len(coords) else 0
        # only include if real edge (length > 0):
        if(coords[c0] != coords[c1]):
            edges.append(LineString([coords[c0], coords[c1]]))
    return edges
#

def checkIfClockwise (edges):
    return sum(list(
        map(lambda x: 
            (list(x.coords)[1][0] - list(x.coords)[0][0])*
            (list(x.coords)[1][1] + list(x.coords)[0][1]), edges)
        )) > 0
#

def isDescending (edgeCoords):
    if(edgeCoords[0][0] == edgeCoords[1][0]):
        return edgeCoords[0][1] > edgeCoords[1][1]
    else: 
        return edgeCoords[0][0] > edgeCoords[1][0]
#

def createLinesFromCoords (coords):
    coords = list(set(coords))
    coords.sort(reverse = isDescending(coords))
    lines = []
    for coordIndex in range(len(coords) - 1):
        newLine = LineString([coords[coordIndex], coords[coordIndex + 1]])
        lines.append(newLine)
    return lines
#

def splitLineWithEdges (line, edges):
    # loop through all edge to evaluate whtether there has been an intersection
    lineCoords = list(line.coords)
    for edge in edges:
        if(line.intersects(edge)):
            # if there has been an intersection cut edge in two parts
            lineCoords += list(edge.intersection(line).coords)
          
    return edges
#

def splitEdgesAlongLine (edges, line):
    oldEdges = edges
    # loop through all edge to evaluate whtether there has been an intersection
    edgesToRemove = []
    edgesToAdd = []
    for edge in edges:
        if(edge.crosses(line)):
            # if there has been an intersection cut edge in two parts
            edgeIntersectionCoords = list(edge.intersection(line).coords)
            edgeCoords = list(set(list(edge.coords) + edgeIntersectionCoords))
            edgesToAdd += createLinesFromCoords(edgeCoords)
            edgesToRemove.append(edge)
        #
    # remove all edges that have been cut in two parts
    for edge in edgesToRemove: edges.remove(edge)
    # add the parts that have been cut in pieces
    edges += edgesToAdd
    
    if(len(edges)<len(oldEdges)-1):
        print("edges removed!")
        print("oldEdges", oldEdges)
        print("edges", edges)
    
    fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
    gpd.GeoSeries(oldEdges).plot(ax=axs, color="red")
    gpd.GeoSeries(edges).plot(ax=axs, color="blue")    
    gpd.GeoSeries(edgesToAdd).plot(ax=axs, color="green")
    gpd.GeoSeries([line]).plot(ax=axs, color="black")    

    return edges
#

def clockwiseAngle (firstAngle, secondAngle):
    return (secondAngle - firstAngle + 540) % 360
#

def plotGeometry (geometryList, axs, color = "", cNr = 0):
        p = gpd.GeoDataFrame(range(len(geometryList)))
        p.columns = ['nr']
        p['geometry'] = geometryList
        cmaps = ['tab20', 'Set3', 'Set1', 'tab20b', 'tab20c']
        if(color):
            if(color not in cmaps):
                p.plot(ax=axs, color=color)
            else:
                p.plot(ax=axs, cmap=color, column='nr')
        else:
            cmap = cmaps[cNr]
            p.plot(ax=axs, cmap=cmap, column='nr')
#
   
for i in range(len(vegShadows.index)):
    # %% section in loop
    row = vegShadows.iloc[i]
    # extract info about polygon
    geom = row['geometry']
    # extract coords of geom. Keep in mind that first point equals last (A,B,C,A)
    coords = list(geom.exterior.coords) # 2do: add coords of interiors.
    
    edges = getEdgesFromCoords(coords)

    # determine if geometry is defined clockwise or counter clockwise
    # clockwise = checkIfClockwise(edges)
    # redefine as clockwise
    # if(not clockwise):
    #     coords.reverse()
    #     # edges = list(map(lambda x: (x[1], x[0]), edges)) #2do adjust to linestring
    #     edges.reverse()
    #     geom.exterior.coords = coords
    #     # possibly interiors might be defined not in the same clock direction as exterior.
    #     for gi in range(len(geom.interiors)):
    #         geom.interiors[gi].coords = list(geom.interiors[gi].coords)[::-1]
    
    
    # calculate angle between sun and earth surface at location and date
    altitude = solar.get_altitude(
        latitude_deg = row['lat4326'], 
        longitude_deg = row['long4326'], 
        when = date
        ) 
    # calculate direction of sunlight at location and date
    azimuth = solar.get_azimuth(
        latitude_deg = row['lat4326'], 
        longitude_deg = row['long4326'], 
        when = date
        )
    
    # length of the shadow
    # 2do how to transform height in meters into shadow offset in degrees / epsg projection
    #   calculate great circle distance for two points in Berlin: center and offset to south so that distance is 1km
    #   do the same for offest to west
    #   divide difference of coordinates by 1000
    #   use this factor to scale in lat and long direction
    #   lat and long might be different
    
    shadowLengthUnitStick = math.tan(altitude * math.pi / 180)
    shadowLength = vegShadows['height'].iloc[i] * 1 * shadowLengthUnitStick
    
    # remember: azimuth gives direction where sun comes from: shadow are projected in the opposite direction
    # if azimuth between 0° and 180° positive (= shadow goes east), else negative
    longShare = - math.sin((azimuth) * math.pi / 180)
    # if azimuth between 90° and 270° negative (= shadow goes south) else negative
    latShare = - math.sin((azimuth - 90) * math.pi / 180)

    xoff = shadowLength * longShare / 1
    yoff = shadowLength * latShare / 1
    print('xoff', xoff, 'yoff', yoff)
    
    geomShifted = shapely.affinity.translate(
        geom = geom, 
        xoff = xoff, 
        yoff = yoff
        )
    coordsShifted = list(geomShifted.exterior.coords)
    edgesShifted = getEdgesFromCoords(coordsShifted)
    # plot polygons
    # plot edges
    fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
    plotGeometry([geom], axs, "#ccccccaa")
    plotGeometry([geomShifted], axs, "#ccccccaa")
    plotGeometry(edges, axs, "blue")
    plotGeometry(edgesShifted, axs, "green")
    fig.suptitle('geom (blue) and geomShifted (green)')

    shifters = []
    intersectionCoords = []
    
    # account for intersections along between geom and geomShifted with Shifters lines
    for coordIndex in range(len(coords)-1):
        shifter = LineString([coords[coordIndex], coordsShifted[coordIndex]])
        shifterCoords = list(shifter.coords)
        if(geom.crosses(shifter)):
            intersectionCoords = list(geom.intersection(shifter).coords)
            shifterCoords += intersectionCoords
            intersectionCoords += intersectionCoords
            # split original edges along shifter 
            edges = splitEdgesAlongLine(edges, shifter)
            
        #
        if(geomShifted.crosses(shifter)):
            intersectionCoords = list(geomShifted.intersection(shifter).coords)
            shifterCoords += intersectionCoords
            intersectionCoords += intersectionCoords
            # split shifted edges along shifter 
            edgesShifted = splitEdgesAlongLine(edgesShifted, shifter)

        # if the shifter has intersections with one of the two polygons
        if(len(shifterCoords) > 2):
            shifterCoords = list(set(shifterCoords))
            shifterCoords.sort(reverse = isDescending(shifterCoords))
            for shifterIndex in range(len(shifterCoords) - 1):
                shifters.append(
                    LineString(
                        [shifterCoords[shifterIndex], 
                         shifterCoords[shifterIndex + 1]]
                        )
                    )
            #
        else:
            shifters.append(shifter)
    #
    

    fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
    plotGeometry(edges, axs, "blue")
    plotGeometry(edgesShifted, axs, "green")
    plotGeometry(shifters, axs, "orange")
    
    # account for intersections between geom and geomShifted
    if(geom.intersects(geomShifted)):
        edgesToRemove = []
        edgesToAdd = []
        for edge in edges:
            if(edge.intersects(geomShifted)):
                # split edge along Edges Shifted
                for edgeShifted in edgesShifted:
                    if(edge.intersects(edgeShifted)):
                        edgesToAdd += splitEdgesAlongLine([edge], edgeShifted)
                edgesToRemove.append(edge)
                # split edgesShifted
                edgesShifted = splitEdgesAlongLine(edgesShifted, edge)
                
            #
        print("remove:", edgesToRemove)
        print("edgesToAdd:", edgesToAdd) # no edges get added...
        for edge in edgesToRemove: edges.remove(edge)
        edges += edgesToAdd
        #
    #
    
    # c0 = [(0,0), (4,0), (4,1), (1,1), (1,3), (0,3)]
    # s0 = Polygon(c0)    
    # s0
    # s1 = Polygon(map(lambda x: (x[0]+3, x[1]+3) , c0))
    # shapely.ops.cascaded_union([s0,s1])    
    # shapely.ops.cascaded_union([s0,s1]).convex_hull
    # l0 = LineString([(0,0),(2,2)])
    # l1 = LineString([(1,1),(3,3)])
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
    
    
    # 2do account for intersections in shifters (=overlaps eg. [(0,0),(2,2)] with [(1,1),(3,3)])
    # might be expensieve
    shiftersToRemove = []
    shiftersToAdd = []
    for shifterIndex in range(len(shifters)):
        shifter = shifters[shifterIndex]
        for otherShifter in (shifters[:shifterIndex] + shifters[shifterIndex+1:]):
            if(shifter.intersects(otherShifter)):
                print("intersection")
                shifter.intersection(otherShifter)
                shiftersToRemove.append(otherShifter)
                shiftersToRemove.append(shifter)
                intersectionCoords = list(shifter.intersection(otherShifter).coords)
                newCoords = intersectionCoords + list(shifter.coords)+ list(otherShifter.coords)
                shiftersToAdd += createLinesFromCoords(newCoords)
    shifterToRemove = list(set(shiftersToRemove))
    for shifter in shifterToRemove: shifters.remove(shifter)
    shifters += shiftersToAdd

    
    fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
    plotGeometry(edges, axs, '', 0)
    plotGeometry(edgesShifted, axs, '', 1)
    plotGeometry(shifters, axs, '', 2)
    fig.suptitle('edges after seperation at intersections')

    
    allEdges = edges + edgesShifted + shifters # somehow this does not include all edges
    fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
    plotGeometry(allEdges, axs, '', 0)
    fig.suptitle('allEdges after seperation at intersections')

    polygonizedEdges = shapely.ops.polygonize_full(allEdges)
    polygonizedEdges = shapely.ops.polygonize_full(edges + edgesShifted + shifters)
    fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
    plotGeometry(polygonizedEdges, axs, '', 0)
    fig.suptitle('allEdges after seperation at intersections')

    allEdgesTwoWayPaths = []
    allEdgeAngles = []
    allPoints = []
        
    for edge in allEdges:
        # are edge unqiue? or only unique considering the direction?
        # i should include inverted directions also
        edgeCoords = list(edge.coords)
        allPoints += edgeCoords
        
        
        # derive angle to north
        xDiff = edgeCoords[0][0] - edgeCoords[1][0]
        yDiff = edgeCoords[1][1] - edgeCoords[0][1]
        angle = 0
        if(xDiff == 0):
            angle = 0 if(yDiff > 0) else 180
        else:
            angle = -math.atan(yDiff / xDiff) * 180 / math.pi + 270 if(xDiff<0) else 90
        
        allEdgeAngles.append(angle)
        allEdgeAngles.append((angle+180) % 360)
        # unfold edges for two way path
        allEdgesTwoWayPaths.append(edgeCoords)
        allEdgesTwoWayPaths.append(edgeCoords[::-1])
        #
    #
    
    allPoints = list(set(allPoints))
    gpd.GeoSeries(
        map(lambda p: shapely.geometry.Point(p),allPoints)
    ).plot(ax=axs, color="000")

    
    shadowPolygons = []
    sunnyPolygons = []
    usedEdgeIndexes = set()
    # function manipulates arrays outside of function scope
    def findNextStep (currentPath):
        candidateArray = list(filter(
            lambda x: 
                # is entry of enumerated array of edges, eg. (0, [(2,2), (4,4)])
                # candidates can not be used yet
                (x[0] not in usedEdgeIndexes) &
                (x[1][0] == allEdgesTwoWayPaths[currentPath[-1]][1]) & 
                # (x[1][0][1] == allEdgesTwoWayPaths[currentPath[-1]][1][1]) &
                # angle in perfectly oposite direction is not valid as it follows same path backwards
                (clockwiseAngle( allEdgeAngles[currentPath[-1]], allEdgeAngles[x[0]] ) > 0), 
                enumerate(allEdgesTwoWayPaths)
            ))
        if(len(candidateArray) == 0):
            print("no candidates")
            return
        #
        
        # idea: we only follow path in clockwise motion
        # then the best way to close polygon is follow the line wich goes 
        # back from which the last edge came from. 
        # the angle thus needs to be maximized but still follow the clockwise definition
        angle = allEdgeAngles[currentPath[-1]]
        bestCandidate = candidateArray[0]
        bestAngle = clockwiseAngle(angle, allEdgeAngles[bestCandidate[0]])
        for candidate in candidateArray[1:]:
            angleWithCandidate = clockwiseAngle(angle, allEdgeAngles[candidate[0]])
            if(angleWithCandidate > bestAngle):
                bestCandidate = candidate
                bestAngle = angleWithCandidate
        
        # add to current path
        currentPath.append(bestCandidate[0])
        # mark edge as used
        usedEdgeIndexes.add(bestCandidate[0])
        # check whether polygon has been closed
        isClosed = allEdgesTwoWayPaths[currentPath[0]][0] == allEdgesTwoWayPaths[currentPath[-1]][1]
        if(not isClosed):
            # find next step
            print("not closed. continue with", bestCandidate[0])
            findNextStep(currentPath)
        else:
            # create and save polygon
            print("closed.")
            newPolygonPoints = []
            for step in currentPath:
                newPolygonPoints.append(
                    allEdgesTwoWayPaths[step][0]
                    )
            newPolygon = Polygon(newPolygonPoints)
            
            # check if new polygon is in shadow
            inShadow = False
            # take a random point from polygon and check if in shadow
            point = shapely.wkt.loads(newPolygon.representative_point().wkt)
            if(geom.contains(point) | geomShifted.contains(point)):
                inShadow = True
            else:
                pointCoords = list(point.coords)[0]
                shifter = LineString(
                    [pointCoords, (pointCoords[0]+xoff, pointCoords[1]+yoff)]
                    )
                if(shifter.intersects(geom)):
                    inShadow = True
                #
            print("inShadow", inShadow)
            newPolygon
            fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
            plotGeometry([newPolygon], axs, ["#ffff0033", "#cccccccc"][inShadow])
            gpd.GeoSeries(map(lambda x: LineString(allEdgesTwoWayPaths[x]),currentPath)).plot(ax=axs, color='greens')
            fig.suptitle('new polygon in sunlight (yellow) or shadow (grey)')

            if(inShadow):
                print("in shadow")
                shadowPolygons.append(newPolygon)
            else:
                print("in sun")
                sunnyPolygons.append(newPolygon)
            #
        #
    #        
    
    for edgeIndex in range(len(allEdgesTwoWayPaths)):
        if(edgeIndex in usedEdgeIndexes):
            continue
        print("try", edgeIndex)
        usedEdgeIndexes.add(edgeIndex)
        #edgeCoords = allEdgesTwoWayPaths[edgeIndex]
        currentPath = [edgeIndex]
        # use nested function to find all paths
        findNextStep(currentPath)
    
    print(shadowPolygons)
    len(shadowPolygons)
    len(sunnyPolygons)
    gpd.GeoSeries(shadowPolygons[:1]).plot(color="#bbbbbb44")
    fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
    gpd.GeoSeries(shadowPolygons).plot(ax=axs, color="#cccccccc")
    gpd.GeoSeries(sunnyPolygons).plot(ax=axs, color="#ffff0033")
    gpd.GeoSeries(edges).plot(ax=axs, color="blue")
    gpd.GeoSeries(edgesShifted).plot(ax=axs, color="green")    
    gpd.GeoSeries(shifters).plot(ax=axs, color="orange")
    gpd.GeoSeries(
        map(lambda p: shapely.geometry.Point(p),allPoints)
    ).plot(ax=axs, color="000")
    fig.suptitle('final shadow and sunlight polygons')

    
    gpd.GeoSeries(shadowPolygons).plot()
    shadow = shapely.ops.cascaded_union(shadowPolygons)
    
    
   # %% section after loop

    # 2do: connect footprint and shadow

    #       alternative approach: 
    #           calculate all intersections with lines
    #           dissolve geometry into these individual shapes
    #               how to do that?
    #               collect all vertices of geom, geomShifted, and Shifters
    #               loop over points: move to point with lowest degree angle to north
    #               for each step save direction aware vertices to usedVertices list
    #               next step shall be lowest angle clockwise relative to last vertice
    #               continue until either:
    #                   original point is reached -> save path as new polygon
    #                   a vertice already is in usedVertices -> discard path
    #               move to next point
    
    # handle edge case: shifters overlap: [(0,0), (2,2)] and [(1,1), (3,3)]
    
    # includes all Points from geom, geomShifted and all intersctions points of edges of 
    # geom, geomShifted and also Shifters
    
    
    #           check one point out of each shape to see if in shadow
    #               shapes coming from shifted footprint dont need to be checked
    #               check by projecting that point back towards the sun and evaluate whether it crosses original footprint
    #           merge all shadow ares to one (multi)-polygon


    # Polygon([[0,0],[2,0], [2,2], [1,2], [1,1],[1,2],  [0,2]]) # invalid
    polygon = Polygon([[1,0],[2,5], [4,4], [4,0]])
    path = LineString([(3,3), (5,5)])
    path.intersects(polygon)
    pathInt = list(path.intersection(polygon).coords)
    comb = list(path.coords)+pathInt
    comb.sort()
    comb
    # geom = polygon
    
    # path = LineString([(1,0), (4,3)])
    # path2 = LineString([(2,1), (3,2)])
    # path.intersects(path2)

    # rotation not necessary to find left, right and bottom
    #       calculate signed difference of points to line
    #       distance parallel to azimuth angle for left and right most
    #       distance to orthogonal line to azimuth angle for bottom point
    
    geomRotated = shapely.affinity.rotate(geom, azimuth)
    coordsRotated = {"xy": list(geomRotated.exterior.coords)}
    coordsRotated["x"] = list(map(lambda x: x[0], coordsRotated["xy"]))
    coordsRotated["y"] = list(map(lambda x: x[1], coordsRotated["xy"]))
    # 2do: check wheter x and y might need to be swapped
    leftIndex = coordsRotated["x"].index(min(coordsRotated["x"]))
    rightIndex = coordsRotated["x"].index(max(coordsRotated["x"]))
    bottomIndex = coordsRotated["y"].index(min(coordsRotated["y"]))
    
    leftIndex < bottomIndex
    rightIndex < bottomIndex
    leftIndex < rightIndex
    plt.plot(coordsRotated["x"], coordsRotated["y"])
    
    # check which way one has to go around.
    if(leftIndex < rightIndex):
        coordsShadow = coords[rightIndex:len(coords)-1] + coords[:leftIndex+1] + coordsShifted[leftIndex:rightIndex+1] + [coords[rightIndex]]
    else:
        coordsShadow = coords[rightIndex:len(coords)-1] + coords[:leftIndex+1] + coordsShifted[leftIndex:rightIndex+1]
       
    shadow = Polygon(coordsShadow)
    shadow
    print(shadow.is_valid)
    print(len(shadow.exterior.coords))
    plt.plot(list(map(lambda x: x[0], coordsShadow[1:12])), list(map(lambda x: x[1], coordsShadow[1:12])) , "r+")
    plt.plot(list(map(lambda x: x[0], shadow.exterior.coords)), list(map(lambda x: x[1], shadow.exterior.coords)))
    plt.plot(list(map(lambda x: x[0], shadow.exterior.coords)), list(map(lambda x: x[1], shadow.exterior.coords)))
    plt.plot(list(map(lambda x: x[0], coords)), list(map(lambda x: x[1], coords)))
    plt.plot(list(map(lambda x: x[0], coordsShifted)), list(map(lambda x: x[1], coordsShifted)))
     
     
    list(map(lambda x: x in coordsShadow, coords))
    list(map(lambda x: x in coords, coordsShadow))
    list(map(lambda x: type(x), coords))
     
    print(len(geom.exterior.coords))
    print(len(geomShifted.exterior.coords))
    
    shadow.difference(geom)
    
    list(map(x in coordsShadow, coords))
    # start at left go over bottom to right#
    # then go to projection at index of right most
    # follow counter clockwise until index of left most
    # close form
    
    # mind the gap: do not index outside bounds...
       
    

        
    # shadow = Polygon([[1,0],[2,5],[4,0]])
    
    connectors = []
    for p in range(len(coords)-1): # - 1 as last entry equals first
        connector =  LineString([
            coords[p], 
            (coords[p][0]+xoff, coords[p][1]+yoff)
            ])
        connectors.append( connector )
        # check whether connector line crosses the polygon
        if(p in [leftIndex, rightIndex, bottomIndex]):
            print("cntn")
            continue
        if(connector.crosses(geom)):
            print("crosses")
            
            shadow = shadow # update shadow
            
        
    print(connectors)
    # are holes in the shadow? footprint+shadow?
    
    # selectNth([1,2,3,4], 2)
    x = [[1,2],[2,3],[3,4]]
    res = map(lambda y: y[0],x)
    print(list(res))
    l = [1,2]
    l.append(3)
    print(l)
    
    vegShadows['geometry'].iloc[i] = shadow


# reproject
vegShadows = vegShadows.to_crs(epsg = crs_code)

fig, axs = plt.subplots(nrows = 1, ncols = 2, figsize = (12,8))
veg.plot(ax=axs[0])
veg.plot(ax=axs[1], cmap = 'Greens', edgecolor = 'black', column = 'Mean_nDOM')
vegShadows.plot(ax=axs[0], color="#cccccccc", edgecolor = 'black')
vegShadows.plot(ax=axs[1], color="#cccccccc", edgecolor = 'black')

# coords = veg.loc[19, 'geometry'].exterior.coords[:]
# coords[0]
# coordList = list(coords[0])
# coordList[0] += 10000
# coordList[1] += 100000
# coords[0] = tuple(coordList)

start=datetime.datetime.now()
X = 100*1000
for x in range(X):
    shadow = Polygon([[1,0],[2,5],[4,0]])

print(datetime.datetime.now()-start, "for ", X, " iterations")

# %% connect shadow projections to footprints

# %% merge and simplify polygons 


# Intersecting Layers intersection = gpd.overlay(layer1, layer2, how = 'intersection')
# Union Layers union = gpd.overlay(layer1, layer2, how = 'union')
# how = [
#   'symmetric_difference'= A|B & !A&B, 
#   'difference' = A & !B, 
#   'dissolve' = merge, 
#       first give union common variable with same value union['common_column'] = 1
#       disolve: union_desolved = union.dissolve(by='commun_column')

# shapely:: object.simplify(tolerance, preserve_topology=True) Returns a simplified representation of the geometric object.

# https://geopandas.org/docs/user_guide/geometric_manipulations.html
# https://geopandas.org/docs/user_guide/data_structures.html
# districts_in_aoi['area'] = districts_in_aoi.area/1000000

# Exporting GeoPandas GeoDataFrames into an ESRI Shapefile
# districts_in_aoi.to_file('districts_within_aoi.shp', driver = "ESRI Shapefile")