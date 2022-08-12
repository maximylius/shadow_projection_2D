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

def orderPresentInArray (el1, el2, array):
    """
    check if el1 and el2 are included in array AND
    that el2 directly follows el1
    including possibility of [el2, ...,el1]
    returns bool
    """
    return len(
        [x for x in array if 
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

def findNextVertex (currentPath, polygonPaths, verticesNetworkMatrix, vertices):
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
        isClosed = nextVertex == currentPath[0]
        if(not isClosed):
            currentPath.append(nextVertex)
            findNextVertex(currentPath, polygonPaths, verticesNetworkMatrix, vertices)
        else:
            # append path to polygonPahts
            polygonPaths.append(currentPath)
    
    # idea: we only follow path in clockwise motion
    # then the best way to close polygon is follow the line wich goes 
    # back from which the last edge came from. 
    # the angle thus needs to be maximized but still follow the clockwise definition
#



# %% project shadow

def projectShadow (geom, height, lat, long, date, method = "manual", showPlots = False):
    """
    return shadow geometry
    takes in original geometry, its height, lat, long, date-time
    method: 
        "difference" (fastest) adds buffers to edges ad substracts these polygons from hull using shapely
        "polygonize" (2nd fastest) uses shapely's 'polygonize' to derive individual polygons
        "manual" (slowest) manually programmed solution to detect individual polygons
    """

    # get list of coords
    coords = list(geom.exterior.coords)
    # get edges for ext
    edges = getEdgesFromCoords(coords)
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
    shadowLengthUnitStick = 1 * math.tan( (90 - altitude) * math.pi / 180)
    shadowLength = height * shadowLengthUnitStick
    
    # remember: azimuth gives direction where sun comes from: 
    # shadow are projected in the opposite direction
    # if azimuth between 0° and 180° positive (= shadow goes east), else negative
    longShare = - math.sin((azimuth) * math.pi / 180)
    # if azimuth between 90° and 270° negative (= shadow goes south) else negative
    latShare = - math.sin((azimuth - 90) * math.pi / 180)
    
    scaleMetersToDegrees = 1 # projection already is in meters
    xoff = shadowLength * longShare / scaleMetersToDegrees
    yoff = shadowLength * latShare / scaleMetersToDegrees
    
    # shift footprint along shadow line (=shifter)
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
    if(showPlots):
        fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
        plotGeometry([geom], axs, "#ccccccaa")
        plotGeometry([geomShifted], axs, "#ccccccaa")
        plotGeometry(edges, axs, "blue")
        plotGeometry(edgesShifted, axs, "green")
        plotGeometry(shifters, axs, "orange")
        fig.suptitle('geom (blue), geomShifted (green) and shifters (orange)')
    
    # divide edges at each intersction with other edges
    if(method != "difference"):
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
                        #
                    #
                    edgesToAdd += addEdges
                    edgesToRemove.append(edge)
                    # split edgesShifted
                    edgesShifted = splitEdgesAlongLine(edgesShifted, edge)
                
                #
            edges = [x for x in edges if not x in edgesToRemove] + edgesToAdd
            #
        #

        # account for intersections in shifters (=overlaps eg. [(0,0),(2,2)] with [(1,1),(3,3)])
        # might be expensieve. do only if offset
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
        
        # if(showPlots):
        #     n=0
        #     fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
        #     plotGeometry(edges, axs, '', 0)
        #     plotGeometry(edgesShifted, axs, '', 1)
        #     plotGeometry(shifters, axs, '', 5)
        #     plotGeometry( list(map(
        #         lambda x: Point(x), 
        #         list(map(lambda x: list(x.coords)[0], shifters)) + list(map(lambda x: list(x.coords)[1], shifters))
        #         )), axs, '#00ffff00', 0)
        #     plotGeometry([edges[n]], axs, 'black', 5)
        #     plotGeometry([Point(edges[n].coords[0]), Point(edges[n].coords[1])], axs, '#ff0000', 5)
        #     fig.suptitle('edges after seperation at intersections (incl shifter overlap)' + str(n))
        #     n += 1 
    #
    
    allEdges = edges + edgesShifted + shifters 

    
    
    precision = 9 # round edge coordinates to prevent precision issues
    polygons = []
    if(method == "difference"):
        # convert edges to MultiLinesString 
        multiLine = shapely.geometry.MultiLineString(allEdges)
        ## Convert Lines to Polygons by applying a tiny buffer
        buf = 0.000001
        multiLine = multiLine.buffer(buf)
        ## Get outer boundary of the lines as a polygon
        boundary = multiLine.convex_hull
        ## Get a difference to generate a multipolygon
        polygons = list(boundary.difference(multiLine))
        # reapply buffer to not have cuts in the shape
        polygons = list(map(lambda x: x.buffer(buf*2), polygons))
        #
    elif(method == "polygonize"):
        allEdgesRounded = list(map(
        lambda x: LineString(list(map(
            lambda x: (round(x[0],precision), round(x[1],precision)),
            list(x.coords)
            ))),
        allEdges
        ))
    
        polygons = list(shapely.ops.polygonize(allEdgesRounded))
        
        # fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
        # plotGeometry(list(polygons), axs, '#22cccc55', 0, True)
        # plotGeometry(allEdgesRounded, axs, '', 5, False)
        # fig.suptitle('allEdges after seperation at intersections')
    elif(method == 'manual'):
        # create matrix of vertex network
        vertices = []
        for edge in allEdges:
            vertices += list(map(
                lambda x: (round(x[0],precision),round(x[1],precision)), 
                edge.coords))
            # vertices += list(edge.coords)
        vertices = list(set(vertices))
        
        edgesOfVertices = [set() for _ in range(len(vertices))]
        for e in range(len(allEdges)): 
            for v1 in range(len(vertices)):
                if(vertices[v1] in list(map(
                        lambda x: (round(x[0],precision),round(x[1],precision)), 
                        allEdges[e].coords))):
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
        
        polygonPaths = []
        for v0 in range(len(vertices)):
            connectedVertices = [x[0] for x in [x for x in enumerate(verticesNetworkMatrix[v0]) if x[1]]]
            for v1 in connectedVertices:
                usedBefore = orderPresentInArray(v0, v1, polygonPaths)
                if(not usedBefore):
                    findNextVertex([v0, v1], polygonPaths, verticesNetworkMatrix, vertices)
        
        
        for polyPath in polygonPaths:
            poly = Polygon(list(map(lambda x: vertices[x], polyPath )))
            polygons.append(poly.buffer(0))

        # fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
        # plotGeometry(list(polygons), axs, '#22cccc55', 0, True)
        # plotGeometry(allEdges, axs, '', 5, False)
        # fig.suptitle('allEdges after seperation at intersections')
    else:
        raise ValueError('unknown method. use manual, polygonize or difference.')
    #
    
    
    shadowPolygons = []
    sunnyPolygons = []
    
    for poly in polygons:
        
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
            shadowPolygons.append(poly)
        else:
            sunnyPolygons.append(poly)
        #

    
    # merge individual shadow polygons to one shadow geometry
    shadow = shapely.ops.cascaded_union(shadowPolygons)
    
    correctlySpecified = (shadow.covers(geom) and shadow.covers(geomShifted)) or (
        shapely.ops.cascaded_union([geom.difference(shadow), geomShifted.difference(shadow)]).area < geom.area / 10000 # allow minor tolerance
        )
    if(not correctlySpecified):
        shadow = shapely.ops.cascaded_union([shadow, geom, geomShifted])
    
    if(showPlots):
        fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))
        gpd.GeoSeries([shadow]).plot(ax=axs, color = ['#ff7777', '#aaffaa'][correctlySpecified])
        if(len(sunnyPolygons)): gpd.GeoSeries(sunnyPolygons).plot(ax=axs, color="#ffff0033")
        gpd.GeoSeries(edges).plot(ax=axs, color="blue")
        gpd.GeoSeries(edgesShifted).plot(ax=axs, color="green")    
        gpd.GeoSeries(shifters).plot(ax=axs, color="orange")
        fig.suptitle('final shadow and sunlight polygons')
    
    return (shadow)
    


# %% loop through

# Month day combinations
filenames = ["Gebaeude_Hoehen_Phase1", "Gebaeude_Hoehen_Phase2"]
dates = [[12,21], [3,21], [6,21], [9,21]]
hours = [18, 9]

for filename in filenames:
    print("open ", filename)
    connection = r'./GeoData/' + filename + '.shp'
    # Importing an ESRI Shapefile
    footprint_df = gpd.read_file(connection)
    print("file read.")
    if("Class_name" in footprint_df.columns):
        footprint_df = footprint_df[
            (footprint_df['Mean_nDOM'] > 0) & 
            (footprint_df['Class_name'] != 'Veg < 2,5 m (einschl. Rasen, Grünland)')
        ]
    
    footprint_df['height'] = footprint_df['Mean_nDOM']

    #reproject to obtain proper WGS84 latitude-longitude projection EPSG:4326, 
    # to obtain correct coordinates to retrieve altitude and azimuth
    footprint_df = footprint_df.to_crs(epsg = 4326)
    footprint_df['centroid'] = footprint_df.centroid
    footprint_df['long4326'] = footprint_df['centroid'].x
    footprint_df['lat4326'] = footprint_df['centroid'].y
    
    # reprojecting to DHDN / Soldner Berlin EPSG:3068 with transformation: 1777
    # DHDN: Deutsches Hauptdreiecksnetz
    # https://epsg.io/3068-15949
    footprint_df = footprint_df.to_crs(epsg = 3068)
    
    # plot using GeoPandas and Matplotlib
    fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (24,16))
    footprint_df.plot(ax=axs, cmap = 'Reds', edgecolor = 'black', column = 'Mean_nDOM')
    fig.suptitle(filename + 'color according to mean height.')


    
    for mmdd in dates:
        for hh in hours:
            shadow_df = footprint_df
            # year, month, day, hour, minute, second, microsecond
            date = datetime(2021, *mmdd, hh, 0, 0, 0, tzinfo=pytz.timezone('Europe/Berlin'))
            altitude = solar.get_altitude(
                latitude_deg = 52.5, 
                longitude_deg = 13.4, 
                when = date
                )
            print("Date:" ,date, "altitude:", altitude)
            if(altitude > 0):
                progress = list(map(lambda x: (len(footprint_df)/100*x)//1, range(1,101)))
                startTime = lastTime =  datetime.now()
                leftOvers = footprint_df.index[i:]
                for i in leftOvers:
                # for i in footprint_df.index:
                    geom = footprint_df.loc[i, 'geometry']#, 
                    height = footprint_df.loc[i, 'height']#,
                    lat = footprint_df.loc[i, 'lat4326']#, 
                    long = footprint_df.loc[i, 'long4326']#
                    method = "difference"#,
                    showPlots = False
                    
                    if(type(geom) != shapely.geometry.polygon.Polygon):
                        continue
                    #
                    
                    (shadow) = projectShadow(
                        geom = geom, 
                        height = height,
                        lat = lat, 
                        long = long, 
                        date = date,
                        method = method,
                        showPlots = showPlots
                        )
                       
                    shadow_df.loc[i, "geometry"] = shadow
                    if(i in progress):
                        print(progress.index(i), "%. Time: ",  datetime.now() - lastTime)
                        lastTime = datetime.now()
                #
                # save as new file
                newFileName = str(filename)+"_"+str(mmdd[0])+"-"+str(mmdd[1])+"_"+str(hh)+".shp"
                print("save to ", newFileName)
                shadow_df.to_file(newFileName)
                
                oneShadowGeom = cascaded_union(shadow_df["geometry"])
                
                # plot using GeoPandas and Matplotlib
                fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (24,16))
                shadow_df.plot(ax=axs, color =  "#ccccccaa")
                footprint_df.plot(ax=axs, cmap = 'Reds', edgecolor = 'black', column = 'Mean_nDOM')
                fig.suptitle(filename + 'color according to mean height.')
                print('Overall time for shadow projection: ', datetime.now() - startTime)
        
# end loop
# print(np.array(perf["manual"]).mean()) # when not silent 0.81s # when silent 0.51s
# print(np.array(perf["polygonize"]).mean()) # when not silent 0.57 # when silent 0.15s
# print(np.array(perf["difference"]).mean()) # when not silent 0. # when silent 0.11s


# %% merge and simplify polygons 

for i in footprint_df.index:
    if(i in progress):
        print(progress.index(i), "%. Time: ",  datetime.now() - lastTime)
        lastTime = datetime.now()
#
                #
                
a = list(map(lambda x: type(x), shadow_df['geometry']))
from itertools import groupby
[len(list(group)) for key, group in groupby(a)]

from shapely.geometry import mapping, Polygon
import fiona

# Here's an example Shapely geometry
df = gpd.GeoSeries([geom])
df.crs = "EPSG:3068"
df = df.to_crs(4326)
poly = geom#Polygon([(0, 0), (0, 1), (1, 1), (0, 0)])

# Define a polygon feature geometry with one attribute
schema = {
    'geometry': 'Polygon',
    'properties': {'id': 'int'},
}

# Write a new Shapefile
with fiona.open('./out/'+newFileName+'ONE'+'.shp', 'w', 'ESRI Shapefile', schema) as c:
    ## If there are multiple geometries, put the "for" loop here
    c.write({
        'geometry': mapping(oneShadowGeom),
        'properties': {'id': newFileName},
    })

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