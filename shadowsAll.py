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
path = 'C:\\Users\\BSE\\Documents\\MaxVonMylius\\GeographyGrowth\\TermPaper\\py'
# path = 'C:\\Users\\mmyli\\Documents\\HumboldtUni\\GeographyGrowth\\Paper\\shadows'
os.chdir(path)

import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import shapely
from shapely.geometry import Polygon, LineString, MultiLineString, Point
from shapely.affinity import translate, scale
from shapely.ops import cascaded_union
from pysolar import solar
from datetime import datetime
import math
import pytz

# %% load shapefile

# Importing an ESRI Shapefile
veg = gpd.read_file(r'./../Sample/sampleVeg.shp')
geb = gpd.read_file(r'./../Sample/sampleGeb.shp')
# veg = gpd.read_file(r'./Sample/sampleVeg.shp')
# geb = gpd.read_file(r'./Sample/sampleGeb.shp')


# check projection
print(veg.crs) # EPSG: 25833
print(geb.crs) # EPSG: 25833

# reprojecting to DHDN / Soldner Berlin EPSG:3068 with transformation: 1777
# DHDN: Deutsches Hauptdreiecksnetz
# https://epsg.io/3068-15949
crs_code = 3068
# project to
veg = veg.to_crs(epsg = crs_code)
geb = geb.to_crs(epsg = crs_code)



# plot using GeoPandas and Matplotlib
fig, axs = plt.subplots(nrows = 2, ncols = 2, figsize = (12,8))
veg.plot(ax=axs[0,0])
veg.plot(ax=axs[0,1], cmap = 'Greens', edgecolor = 'black', column = 'Mean_nDOM')
geb.plot(ax=axs[1,0], cmap = 'Reds', edgecolor = 'black', column = 'Mean_nDOM')
geb.plot(ax=axs[1,1], cmap = 'Reds', edgecolor = 'black', column = 'Max_nDOM')
veg.plot(ax=axs[0,0], color = "#442288", edgecolor = "pink")


geb
geb.iloc[20, list(geb.columns).index('geometry')]

# %% set up sunlight angles


# year, month, day, hour, minute, second, microsecond
# date = datetime(2007, 2, 16, 12, 13, 1, 130320, tzinfo=pytz.timezone('Europe/Berlin'))


# latitude_southMostBerlin = 52.388447
# longitude_westMostBerlin = 13.131997
# latitude_northMostBerlin = 52.675258
# longitude_northMostBerlin = 13.479466
# print(solar.get_azimuth(latitude_southMostBerlin, longitude_westMostBerlin, date))
# print(solar.get_azimuth(latitude_northMostBerlin, longitude_northMostBerlin, date))

# print(solar.get_altitude(latitude_southMostBerlin, longitude_westMostBerlin, date))
# print(solar.get_altitude(latitude_northMostBerlin, longitude_northMostBerlin, date))


# latitude = 42.206
# longitude = -71.382



# print(date)
# # calculate angle between sun and earth surface at location and date
# altitude = solar.get_altitude(latitude, longitude, date) 
# # calculate direction of sunlight at location and date
# azimuth = solar.get_azimuth(latitude, longitude, date)
# print(altitude)
# print(azimuth)

# if azimuth between 0° and 180° positive (=going east), else negative
# x = math.sin((azimuth) * math.pi / 180)
# if azimuth between 270° and 90° positive (=going north) else negative
# y = math.sin((azimuth - 270) * math.pi / 180)

# shadowLength = math.tan(altitude * math.pi / 180)]


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
vegShadows
# project back to berlin projection
vegShadows = vegShadows.to_crs(epsg = crs_code)


vegShadows['geomShifted'] = vegShadows['geometry']
vegShadows['shadowProjected'] = vegShadows['geometry']
vegShadows['shadow'] = vegShadows['geometry']


date = datetime(2007, 4, 16, 16, 13, 1, 130320, tzinfo=pytz.timezone('Europe/Berlin'))

# %% helper functions outside scope
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

def splitEdgesAlongLine (edges, line, plot = False):
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
    
    if(plot):
        fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
        gpd.GeoSeries(oldEdges).plot(ax=axs, color="red")
        gpd.GeoSeries(edges).plot(ax=axs, color="blue")    
        gpd.GeoSeries(edgesToAdd).plot(ax=axs, color="green")
        gpd.GeoSeries([line]).plot(ax=axs, color="black")
        fig.suptitle('Spliting (green) edge at (black) line')    

    return edges
#
def angleToNorth (point1 = (0,0), point2 = (0,0)):
        """
        point1 and point2 each must be a tuple of point coords eg. (0,1)
        return angle to north ((0,0) to (0,1))
        """
        xDiff = point2[0] - point1[0]
        yDiff = point2[1] - point1[1]
        
        if(xDiff == 0):
            angle = 0 if(yDiff > 0) else 180
        else:
            angle = -math.atan(yDiff / xDiff) * 180 / math.pi + (270 if(xDiff<0) else 90)
        
        return angle
#    

def clockwiseAngle (firstAngle, secondAngle):
    """
    returns the outer angle in clockwise direction
    thus 180° if angles are equal
    """
    return (secondAngle - firstAngle + 540) % 360
#

def plotGeometry (geometryList, axs = None, color = "", cNr = None, label = False):
    
    # """
    # plot a list of geometries. input:
    #     geometryList = list of geometries
    #     axs = plot index, if not specified it directly plots
    #     color = if specified it determines the color of geometry
    #     cNr = (if no color) species color map ['tab20', 'Set3', 'Set1', 'tab20b', 'tab20c', 'gist_ncar'][cNr].  
    # """
    
    if( type(geometryList) != list or (len(geometryList) == 0) ):
        print('can not plot geometry: no valid list of geometries provided.')
        return
    
    p = gpd.GeoDataFrame(range(len(geometryList)))
    p.columns = ['nr']
    p['geometry'] = geometryList
    cmaps = ['tab20', 'Set3', 'Set1', 'tab20b', 'tab20c', 'gist_ncar']
    if(color):
        if(color not in cmaps):
            p.plot(ax = axs, color = color, edgecolor = '#000')
        else:
            p.plot(ax = axs, cmap = color, column = 'nr')
    else:
        cmap = cmaps[cNr]
        p.plot(ax = axs, cmap = cmap, column = 'nr')
    if(axs and label):
        p.apply(
            lambda x: axs.annotate(text=x.nr, xy=x.geometry.centroid.coords[0], ha='center'), axis=1
        );
#

def projectShadow (geom, height, lat, long, date):
    
    # geom = vegShadows.loc[i, 'geometry']#, 
    # get list of coords
    coords = list(geom.exterior.coords)
    # coords = coords[27:36]
    # get edges for ext
    edges = getEdgesFromCoords(coords)
    # gpd.GeoSeries(edges).plot()
    # list(shapely.ops.polygonize(edges))
    # gpd.GeoSeries(list(shapely.ops.polygonize(edges))).plot()
    geom = shapely.ops.cascaded_union(list(shapely.ops.polygonize(edges)))
    # add edges and coordinates of possible interiors (holes in shape)
    for interior in list(geom.interiors):
        interiorCoords = list(interior.coords)
        interiorEdges = getEdgesFromCoords(interiorCoords)
        coords += interiorCoords
        edges += interiorEdges
    
    # calculate angle between sun and earth surface at location and date
    altitude = solar.get_altitude(
        latitude_deg = lat, 
        longitude_deg = long, 
        when = date
        ) 
    # calculate direction of sunlight at location and date
    azimuth = solar.get_azimuth(
        latitude_deg = lat, 
        longitude_deg = long, 
        when = date
        )
    # length of the shadow
    shadowLengthUnitStick = math.tan(altitude * math.pi / 180)
    shadowLength = height * 1 * shadowLengthUnitStick
    
    # remember: azimuth gives direction where sun comes from: 
    # shadow are projected in the opposite direction
    # if azimuth between 0° and 180° positive (= shadow goes east), else negative
    longShare = - math.sin((azimuth) * math.pi / 180)
    # if azimuth between 90° and 270° negative (= shadow goes south) else negative
    latShare = - math.sin((azimuth - 90) * math.pi / 180)
    
    scaleMetersToDegrees = 1 # projection already is in meters
    xoff = shadowLength * longShare / scaleMetersToDegrees
    yoff = shadowLength * latShare / scaleMetersToDegrees
    print('y/x', yoff/xoff, 'xoff', xoff, 'yoff', yoff)
    
    geomShifted = translate(
        geom = geom, 
        xoff = xoff, 
        yoff = yoff
        )
    
    coordsShifted = list(geomShifted.exterior.coords)
    edgesShifted = getEdgesFromCoords(coordsShifted)
    # add edges and coordinates of possible interiors (holes in shape)
    for interior in list(geomShifted.interiors):
        interiorCoords = list(interior.coords)
        interiorEdges = getEdgesFromCoords(interiorCoords)
        coordsShifted += interiorCoords
        edgesShifted += interiorEdges
        
    # get shifters (connecting point of original geometry to respective point in projection)
    shifters = []    
    for coordIndex in range(len(coords)-(coords[0]==coords[-1])):
        shifter = LineString([coords[coordIndex], coordsShifted[coordIndex]])
        shifters.append(shifter)

    # plot polygons, edges and shifters
    fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
    plotGeometry([geom], axs, "#ccccccaa")
    plotGeometry([geomShifted], axs, "#ccccccaa")
    plotGeometry(edges, axs, "blue")
    plotGeometry(edgesShifted, axs, "green")
    plotGeometry(shifters, axs, "orange")
    fig.suptitle('geom (blue), geomShifted (green) and shifters (orange)')
        
    # sect intersection of shifters with geoms
    intersectionCoords = []
    shiftersToRemove = []
    shiftersToAdd = []
    # account for intersections along between geom and geomShifted with Shifters lines
    for coordIndex in range(len(coords)-1):
        
        shifter = shifters[coordIndex]
        shifterCoords = list(shifter.coords)
        
        if(geom.crosses(shifter)):
            # print('crosses geom')
            # there is a mistake here. intersection can be point only
            intersections = geom.intersection(shifter)
            # print('intersections.type' , intersections.type)
            if(intersections.type in ['LineString', 'Point']):
                intersections = [intersections]
            for intersection in intersections: 
                # if(intersection.type == 'Point'):
                #     print('point only')
                #     continue
                intersectionCoords = list(intersection.coords)
                shifterCoords += intersectionCoords
                intersectionCoords += intersectionCoords
                # split original edges along shifter 
                edges = splitEdgesAlongLine(edges, shifter)
            
        #
        if(geomShifted.crosses(shifter)):
            # print('crosses geomShifted')
            intersections = geomShifted.intersection(shifter)
            # print('intersections.type' , intersections.type)
            if(intersections.type in ['LineString', 'Point']):
                intersections = [intersections]
            for intersection in intersections: 
                # if(intersection.type == 'Point'):
                #     print('point only')
                #     continue
                intersectionCoords = list(intersection.coords)
                shifterCoords += intersectionCoords
                intersectionCoords += intersectionCoords
                # split shifted edges along shifter 
                edgesShifted = splitEdgesAlongLine(edgesShifted, shifter)

        # if the shifter has intersections with one of the two polygons        
        shifterCoords = list(set(shifterCoords))
        if(len(shifterCoords) > 2):
            shiftersToRemove.append(shifter)
            shifterCoords.sort(reverse = isDescending(shifterCoords))
            for shifterIndex in range(len(shifterCoords) - 1):
                # print("add 1 shifter")
                shiftersToAdd.append(
                    LineString(
                        [shifterCoords[shifterIndex], 
                         shifterCoords[shifterIndex + 1]]
                        )
                    )
            #
        #
    #
    shifters = [x for x in shifters if x not in shiftersToRemove] + shiftersToAdd

    
    #sec intersection geom & geomShifted
    
    # account for intersections between geom and geomShifted
    if(geom.intersects(geomShifted) and (not geom.touches(geomShifted))):
        edgesToRemove = []
        edgesToAdd = []
        for edge in edges:
            if(edge.crosses(geomShifted)):
                # print("edge crossses---")
                addEdges = []
                # split edge along Edges Shifted
                for edgeShifted in edgesShifted:
                    # if(edge.intersects(edgeShifted)): print(edge.crosses(edgeShifted), edge.intersects(edgeShifted), not edge.touches(edgeShifted))
                    if(edge.intersects(edgeShifted) and (not edge.touches(edgeShifted))):
                        addEdges = splitEdgesAlongLine([edge], edgeShifted)
                
                fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
                plotGeometry([geomShifted], axs, "#cccccc33",0)
                plotGeometry(edges, axs, "blue",1)
                plotGeometry([edge], axs, "red",1)
                plotGeometry(addEdges, axs, "green",1)
                plotGeometry([x for x in splitEdgesAlongLine(edgesShifted, edge) if not x in edgesShifted], axs, "orange",1)

                edgesToAdd += addEdges
                edgesToRemove.append(edge)
                # split edgesShifted
                edgesShifted = splitEdgesAlongLine(edgesShifted, edge)
            
            #
        # print("remove:", edgesToRemove) # removes to many
        # print("edgesToAdd:", edgesToAdd) # not enough edges get added...
        edges = [x for x in edges if not x in edgesToRemove] + edgesToAdd
        #
    #
    
    # fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
    # plotGeometry(edges, axs, "",1)
    # plotGeometry(edgesShifted, axs, "",1)
    # plotGeometry(shifters, axs, "", 2)
    # fig.suptitle('edges after seperation at intersections (incl geom overlap)')
 
    # 2do account for intersections in shifters (=overlaps eg. [(0,0),(2,2)] with [(1,1),(3,3)])
    # might be expensieve
    shiftersToRemove = []
    shiftersToAdd = []
    for shifterIndex in range(len(shifters)):
        shifter = shifters[shifterIndex]
        for otherShifter in (shifters[:shifterIndex] + shifters[shifterIndex+1:]):
            if(shifter.intersects(otherShifter) and (not shifter.touches(otherShifter))):
                # print("intersection")
                shifter.intersection(otherShifter)
                shiftersToRemove.append(otherShifter)
                shiftersToRemove.append(shifter)
                intersectionCoords = list(shifter.intersection(otherShifter).coords)
                newCoords = intersectionCoords + list(shifter.coords)+ list(otherShifter.coords)
                shiftersToAdd += createLinesFromCoords(newCoords)
    shifters = [x for x in shifters if not x in shiftersToRemove] + shiftersToAdd

    n=0
    fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
    plotGeometry(edges, axs, '', 0)
    plotGeometry(edgesShifted, axs, '', 1)
    plotGeometry(shifters, axs, '', 5)
    plotGeometry( list(map(
        lambda x: Point(x), 
        list(map(lambda x: list(x.coords)[0], shifters)) + list(map(lambda x: list(x.coords)[1], shifters))
        )), axs, '#00ffff00', 0)
    plotGeometry([edges[n]], axs, 'black', 5)
    plotGeometry([Point(edges[n].coords[0]), Point(edges[n].coords[1])], axs, '#ff0000', 5)
    fig.suptitle('edges after seperation at intersections (incl shifter overlap)' + str(n))
    n += 1 
    
    # problemPoints = [
    #    edges[12].coords[0],
    #    edges[12].coords[1],
    #    edges[14].coords[0],
    #    edges[14].coords[1]
    #    ]
    # list(map(
    #     lambda x: x in coords or x in coordsShifted,
    #     [edges[10].coords[0], edges[10].coords[1]]
    #     #problemPoints
    #         ))
    
    allEdges = edges + edgesShifted + shifters 
    scaledEdges = list(map(lambda x: scale(x, 2, 2), allEdges))
    # 2do scale them up until bounding box
    # result, dangles, cuts, invalids = shapely.ops.polygonize_full(scaledEdges)    
    # gpd.GeoSeries([result, dangles, cuts, invalids]).plot()
    # result.area
    # dangles.area
    # cuts.area
    # invalids.area
    # p = list(shapely.ops.polygonize([
    #             LineString([(0,0),(1,1)]),
    #             LineString([(0,0),(1,0)]),
    #             LineString([(0,0),(0,1)]),        
    #             LineString([(1,0),(1,1)]),
    #             LineString([(0,1),(1,1)]),
    #             LineString([(-1,2),(0,-1)])
    #                ]))
    # p[0]
    # p[1]
    # allPolygons = list(shapely.ops.polygonize_full(allEdges))
    allPolygons = list(shapely.ops.polygonize(allEdges))
    allPolygons = list(shapely.ops.polygonize(scaledEdges))
    allPolygons = list(shapely.ops.polygonize(MultiLineString(allEdges)))
    # shapely.ops.polygonize_full(allEdges)[0]
    # shapely.ops.polygonize_full(allEdges)[3]
    # list(shapely.ops.polygonize_full(allEdges)[0])
    # list(shapely.ops.polygonize_full(scaledEdges)[3])
    
    fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
    plotGeometry(list(allPolygons), axs, '#22cccc55', 0)
    plotGeometry(allEdges, axs, '', 5, True)
    fig.suptitle('allEdges after seperation at intersections')
    
    
  


    shadowPolygons = []
    sunnyPolygons = []
    
    for poly in allPolygons:
        inShadow = False
        # take a random point from polygon and check if in shadow
        point = shapely.wkt.loads(poly .representative_point().wkt)
        if(geom.contains(point) | geomShifted.contains(point)):
            # print("in original")
            inShadow = True
        else:
            # print("not in original")
            pointCoords = list(point.coords)[0]
            shifter = LineString([
                pointCoords, 
                (pointCoords[0]-xoff, pointCoords[1]-yoff)
                ])
            if(shifter.intersects(geom)):
                inShadow = True
            #
        # print("inShadow", inShadow)
        # fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
        # plotGeometry([poly ], axs, ["#ffff0033", "#cccccccc"][inShadow])
        # fig.suptitle('new polygon in sunlight (yellow) or shadow (grey)')

        if(inShadow):
            # print("in shadow")
            shadowPolygons.append(poly)
        else:
            # print("in sun")
            sunnyPolygons.append(poly)
        #
    # fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
    # plotGeometry([geom], axs, "blue")
    # plotGeometry([geomShifted], axs, "green")
    # plotGeometry(shadowPolygons, axs, "#ccccccdd")
    # plotGeometry(sunnyPolygons, axs, "#ffff0077")
    # fig.suptitle('sunlight (yellow) and shadow (grey)')
    
    
    # create matrix of vertex network
    vertices = []
    for edge in allEdges:
        vertices += list(map(lambda x: (round(x[0],10),round(x[1],10)), edge.coords))
        # vertices += list(edge.coords)
    vertices = list(set(vertices))
    
    edgesOfVertices = [set() for _ in range(len(vertices))]
    for e in range(len(allEdges)): 
        for v1 in range(len(vertices)):
            if(vertices[v1] in list(map(lambda x: (round(x[0],10),round(x[1],10)), allEdges[e].coords))):
                edgesOfVertices[v1].add(e)
    
    verticesNetworkMatrix = np.zeros((len(vertices), len(vertices))) # fill with 0
    for v1 in range(len(vertices)):
        for v2 in range(len(vertices)):
            if(not v2 >= v1 and edgesOfVertices[v1].intersection(edgesOfVertices[v2])):
                verticesNetworkMatrix[v1, v2] = 1 # lower triangle matrix
                verticesNetworkMatrix[v2, v1] = 1 # upper trianle matrix

    # the number of polygons a vertex is included in equals
    # the number of edges connected to the point OR this number minus 1
    # does the vertice Network Matrix tell the number of polygons? no.
    
    # work the way through the matrix
    # how to know whether a point can be inlcuded in another polygon?
    # for each step the angle needs to be calculated
    # check wheter path has been taken already:
    # by checking of the in the vertices path list,
    #   the vertices is included
    #   and if its directly followed by the other edge 
    
    def orderPresentInArray (el1, el2, array):
            """
            check if el1 and el2 are included in array AND
            that el2 directly follows el1
            including possibility of [el2, ...,el1]
            returns bool
            """
            return len(
                [x for x in polygonPaths if 
                 el1 in x and 
                 el2 in x and 
                 (
                     x.index(el1) == x.index(el2) - 1 or
                     (
                         x.index(el1) == len(x) - 1 and 
                         x.index(el2) == 0)
                     )
                 ]
                ) > 0
    #
    fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
    plotGeometry(allEdges, axs, "#eee")
    plotGeometry(
        list(map(lambda x: Point(x), vertices)), axs, "#0000", 0,  True
        )
    plotGeometry(
        list(map(lambda x: Point(x), np.array(vertices)[
            [5]##polygonPaths[15] # [5, 55, 9, 39, 16, 43, 26, 63, 55]
            ] )), axs, "red", 0,  False
        )
    

    def findNextVertex (currentPath, polygonPaths):
        # problem: outline path is still included # idea: could be solved if always going the wrong outside way
        """
        currentPath needs to be at least of length 2
        currentPath entries are indices referring to a vertex
        """
        # all other vertices that are connected to this vertex
        lastVertex = currentPath[-2]
        currentVertex = currentPath[-1]
        connectedVertices = [x[0] for x in [x for x in enumerate(verticesNetworkMatrix[currentVertex]) if x[1]]]
        
        angleLastEdge = angleToNorth(vertices[lastVertex], vertices[currentVertex])
        #
        angles = np.array(list(map(
            lambda x: clockwiseAngle(
                           firstAngle = angleLastEdge, 
                           secondAngle= angleToNorth(vertices[currentPath[-1]], vertices[x])
                           ), connectedVertices)))
        nextVertex = connectedVertices[angles.argmax()]
        # check if nextVertex already followed nextVertex in a polyognPaths
        
        usedBefore = orderPresentInArray(currentVertex, nextVertex, polygonPaths)
        
        if(not usedBefore):
            # check if path is closed
            print(nextVertex , currentPath[0])
            isClosed = nextVertex == currentPath[0]
            if(not isClosed):
                currentPath.append(nextVertex)
                findNextVertex(currentPath, polygonPaths)
                print("nest deeper, path", currentPath)
            else:
                # append path to polygonPahts
                polygonPaths.append(currentPath)
        
        # idea: we only follow path in clockwise motion
        # then the best way to close polygon is follow the line wich goes 
        # back from which the last edge came from. 
        # the angle thus needs to be maximized but still follow the clockwise definition
    #
    
    polygonPaths = []
    for v0 in range(len(vertices)):
        connectedVertices = [x[0] for x in [x for x in enumerate(verticesNetworkMatrix[v0]) if x[1]]]
        for v1 in connectedVertices:
            usedBefore = orderPresentInArray(v0, v1, polygonPaths)
            if(not usedBefore):
                findNextVertex([v0, v1], polygonPaths)
    
    
    shadowPolygons = []
    sunnyPolygons = []
    polygons = []
    for polyPath in polygonPaths:
        poly = Polygon(list(map(lambda x: vertices[x], polyPath )))
        polygons.append(poly)
        
        inShadow = False
        # take a random point from polygon and check if in shadow
        point = shapely.wkt.loads(poly .representative_point().wkt)
        if(geom.contains(point) | geomShifted.contains(point)):
            # print("in original")
            inShadow = True
        else:
            # print("not in original")
            pointCoords = list(point.coords)[0]
            shifter = LineString([
                pointCoords, 
                (pointCoords[0]-xoff, pointCoords[1]-yoff)
                ])
            if(shifter.intersects(geom)):
                inShadow = True
            #
        # print("inShadow", inShadow)
        # fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
        # plotGeometry([poly ], axs, ["#ffff0033", "#cccccccc"][inShadow])
        # fig.suptitle('new polygon in sunlight (yellow) or shadow (grey)')

        if(inShadow):
            # print("in shadow")
            shadowPolygons.append(poly)
        else:
            # print("in sun")
            sunnyPolygons.append(poly)
        #
    
    plotGeometry(polygons, None, "#55555511", 5, True)
    # 6 
    # 21, 11, 12, 14, 15
    polygonPaths[6]

    # merge individual shadow polygons to one shadow geometry
    shadow = shapely.ops.cascaded_union(shadowPolygons)
    
    correctlySpecified = (shadow.covers(geom) and shadow.covers(geomShifted)) or (
        shapely.ops.cascaded_union([geom.difference(shadow), geomShifted.difference(shadow)]).area < geom.area / 10000 # allow minor tolerance
        )
    if(not correctlySpecified):
        shadow = shapely.ops.cascaded_union([shadow, geom, geomShifted])
    
    # fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
    # gpd.GeoSeries([shadow]).plot(ax=axs, color = ['#ff7777', '#aaffaa'][correctlySpecified])
    # gpd.GeoSeries(sunnyPolygons).plot(ax=axs, color="#ffff0033")
    # gpd.GeoSeries(edges).plot(ax=axs, color="blue")
    # gpd.GeoSeries(edgesShifted).plot(ax=axs, color="green")    
    # gpd.GeoSeries(shifters).plot(ax=axs, color="orange")
    # fig.suptitle('final shadow and sunlight polygons')
    
    return (shadow, geomShifted)
    
# %% loop through
j = 10


j += 1
i = vegShadows.index[j]
vegShadows.loc[vegShadows.index[j], 'geometry']

j -= 1


for i in vegShadows.index:
    startTime = datetime.now()
    geom = vegShadows.loc[i, 'geometry']#, 
    height = vegShadows.loc[i, 'height']#,
    lat = vegShadows.loc[i, 'lat4326']#, 
    long = vegShadows.loc[i, 'long4326']#
        
    (shadow, geomShifted) = projectShadow(
        geom = geom, 
        height = height,
        lat = lat, 
        long = long, 
        date = date
        )
       
    vegShadows.loc[i, "geomShifted"] = geomShifted
    vegShadows.loc[i, "shadow"] = shadow
    
    
    print('This iteration took ', datetime.now()-startTime)
# end loop
    
# %% section after loop
fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
plotGeometry(list(vegShadows['shadow']), axs, "#cccccccc")
plotGeometry(list(vegShadows['geometry']), axs, "green")
# plotGeometry(list(veg['geometry']), axs, "green")



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

    # rotation not necessary to find left, right and bottom
    #       calculate signed difference of points to line
    #       distance parallel to azimuth angle for left and right most
    #       distance to orthogonal line to azimuth angle for bottom point
    

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