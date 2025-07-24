import sys
import numpy as np
import pyvista as pv
from scipy.spatial import ConvexHull, Delaunay


def calc_circumcenter_circumsphere_tetrahedron_2(tet):
    '''An alternative implementation based on http://mathworld.wolfram.com/Circumsphere.html because of issues with the initial implementation from the Berkeley page.'''
    p0 = tet[0]
    p1 = tet[1]
    p2 = tet[2]
    p3 = tet[3]

    p0squared = p0[0]*p0[0] + p0[1]*p0[1] + p0[2]*p0[2] 
    p1squared = p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2] 
    p2squared = p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2] 
    p3squared = p3[0]*p3[0] + p3[1]*p3[1] + p3[2]*p3[2] 

    D_x = np.linalg.det(np.array([[p0squared, p0[1], p0[2], 1], 
                                    [p1squared, p1[1], p1[2], 1], 
                                    [p2squared, p2[1], p2[2], 1], 
                                    [p3squared, p3[1], p3[2], 1]]))


    D_y = - np.linalg.det(np.array([[p0squared, p0[0], p0[2], 1], 
                                    [p1squared, p1[0], p1[2], 1], 
                                    [p2squared, p2[0], p2[2], 1], 
                                    [p3squared, p3[0], p3[2], 1]]))

    D_z = np.linalg.det(np.array([[p0squared, p0[0], p0[1], 1], 
                                    [p1squared, p1[0], p1[1], 1], 
                                    [p2squared, p2[0], p2[1], 1], 
                                    [p3squared, p3[0], p3[1], 1]]))

    a = np.linalg.det(np.array([[p0[0], p0[1], p0[2], 1],
                               [p1[0], p1[1], p1[2], 1],
                               [p2[0], p2[1], p2[2], 1],
                               [p3[0], p3[1], p3[2], 1]]))

    denominator = 2. * a
    x_0 = D_x / denominator
    y_0 = D_y / denominator
    z_0 = D_z / denominator
    circumcenter = np.array([x_0, y_0, z_0])
    return circumcenter

class Vertex:

    def __init__(self, p):
        self.position = p
        self.tets = []
        self.faces = []
        self.edges = []
        self.isHullVertex = 0


def findPoint(point, verticesList):
    res = -1
    for i in range(0, len(verticesList)):
        vpos = verticesList[i].position
        if abs(point[0] - vpos[0]) < 1e-9 and abs(point[1] - vpos[1]) < 1e-9 and abs(point[2] - vpos[2]) < 1e-9 : 
            res = i 
            break 
    return res


class Edge:
    

    def __init__(self, edgeObject):
        self.points = edgeObject.points
        self.center = edgeObject.center
        self.vertices = []
        self.faces = []
        self.tets = []
        self.isHullEdge = 0


    def equal(self, anotherEdge):
        return (np.array_equiv(self.start, anotherEdge.start) and np.array_equiv(self.end, anotherEdge.end)) or ( np.array_equiv(self.start, anotherEdge.end) and np.array_equiv(self.end, anotherEdge.start))

def findEdge(edgeObject, edgesList):
    res = -1
    center = edgeObject.center
    for i in range(0, len(edgesList)):
        ecenter = edgesList[i].center
        if abs(center[0] - ecenter[0]) < 1e-9 and abs(center[1] - ecenter[1]) < 1e-9 and abs(center[2] - ecenter[2]) < 1e-9: 
            res = i 
            break 
    return res

class Face:

    def __init__(self, faceObject):
        self.points = faceObject.points
        self.center = faceObject.center

        self.edges = []
        self.tets = []

        self.isHullFace = 0


def findFace(faceObject, facesList):
    res = -1
    center = faceObject.center
    for i in range(0, len(facesList)):
        fcenter = facesList[i].center
        if abs(center[0] - fcenter[0]) < 1e-9 and abs(center[1] - fcenter[1]) < 1e-9 and abs(center[2] - fcenter[2]) < 1e-9:
            res = i 
            break 
    return res


class Tetrahedro:



    def __init__(self, p):
        self.points = p
        self.faces = []

        self.edges = []
        self.circumcenter = -1
        self.isHullTet = 0
        self.hullFace = -1

    def getCircumcenter(self):
        return calc_circumcenter_circumsphere_tetrahedron_2(self.points)
    
    def setOuterFaceCenter(self, p, j):
        self.outerFaceCenter = p
        self.outerFaceIndex = j

def existingEdge(edgeArray, edge):
    for i in range(0, len(edgeArray)):
        if edgeArray[i].equal(edge):
            return (True, i)
        
    return (False, -1)






mesh = pv.read(sys.argv[1])

bordercells = []

newborderVertices = []






tetrahedra = []
tetFaces = []
tetEdges = []
tetVertices = []
tetVerticesPositions = []

voronoiVertices = []

countOFnewFaces = 0

print("first pass trough cells")

for i in range(0, mesh.n_cells):

    cell = mesh.get_cell(i)

    tet = Tetrahedro(cell.points)

    cell_faces = mesh.get_cell(i).faces
    for face in cell_faces:
        f0 = Face(face)

        newFace = 0
        index = findFace(f0, tetFaces)
        if index == -1:
            countOFnewFaces+=1
            newFace = 1
            index = len(tetFaces) #apendiar al final
        else:
            f0 = tetFaces[index]

        edgeObjects = []

        faceEdges = face.edges
        for edge in faceEdges:
            e0 = Edge(edge)

            newEdge = 0
            edgeIndex = findEdge(e0, tetEdges)
            if edgeIndex == -1:
                newEdge = 1
                edgeIndex = len(tetEdges)
            else:
                e0 = tetEdges[edgeIndex]


            vertexObjects = []

            edgeVertices = edge.points
            for vertex in edgeVertices:
                p0 = Vertex(vertex) #vertex is array of 3 points
                point0Index = findPoint(vertex, tetVertices)
                if point0Index == -1:
                    point0Index = len(tetVertices)
                    p0.tets.append(i)
                    p0.faces.append(index)
                    p0.edges.append(edgeIndex)
                    tetVertices.append(p0)
                    tetVerticesPositions.append(np.array(vertex))
                else:
                    p0 = tetVertices[point0Index]
                    if not (i in p0.tets): p0.tets.append(i)
                    if not(index in p0.faces): p0.faces.append(index)
                    if not(edgeIndex in p0.edges): p0.edges.append(edgeIndex)

                vertexObjects.append(point0Index)
            
            for vertex in vertexObjects:
                if not (vertex in e0.vertices): e0.vertices.append(vertex) #add vertices to edge

            if not (index in e0.faces) : e0.faces.append(index)     # add face to edge
            if not (i in e0.tets) : e0.tets.append(i)               #adds tetrahedra

            if newEdge: tetEdges.append(e0)

            edgeObjects.append(edgeIndex)

        for edge in edgeObjects:
            if not (edge in f0.edges): f0.edges.append(edge) #add edges to face

        if not (i in f0.tets): f0.tets.append(i)        # add tet to face

        if newFace: tetFaces.append(f0)             # se agrega a la lista de faces si es que no existe

        tet.faces.append(index)


    #obtiene circumcentro
    tetCircum = tet.getCircumcenter()
    # se guarda el indice que tendrá el vertice
    tet.circumcenter = (len(voronoiVertices), tetCircum)
    # se guarda el actual vertice
    voronoiVertices.append(tetCircum)

    tetrahedra.append(tet)

surfaceFacesObjects = []

pl = pv.Plotter()

# marking surface tets and faces
for i in range(0, mesh.n_cells):                            # por cada tetraedra
    neigh = mesh.cell_neighbors(i, "faces")                 # obtengo los vecinos
    if len(neigh) < 4:      # border tet                    # si tiene menos de 4 vecinos es hull tet, me consta que esto esta bien

        tetrahedra[i].isHullTet = 1

        if len(neigh) == 2:
            #print("aaa")
            #for face in tetrahedra[i].faces: tetFaces[face].isHullFace = 1
            tetFaces[tetrahedra[i].faces[0]].isHullFace = 1 # seteo primera cara como hull no más a ver que pasa
            tetFaces[tetrahedra[i].faces[1]].isHullFace = 1
            tetrahedra[i].hullFace = (tetrahedra[i].faces[0], tetrahedra[i].faces[1])
            continue

        tetHasHullFace = 0
                

        for face in tetrahedra[i].faces:                    # para cada face del tetraedro
            found = 0                                       # digo que no he visto la cara

            for neighbor in neigh:  #neighbor es un indice  # para cada tetraedro vecino
                
                neighFaces = tetrahedra[neighbor].faces     # obtengo sus faces
                #print("comparing: ", tetrahedra[i].faces, " con ", neighFaces)
                if face in neighFaces:                      # reviso si mi cara actual tambien es parte de las caras del vecino
                    #print("found")
                    found = 1
                    break
            
            if not found:
                #print("not found")
                tetHasHullFace = 1
                tetFaces[face].isHullFace = 1
                tetrahedra[i].hullFace = face

                for edge in tetFaces[face].edges: 
                    tetEdges[edge].isHullEdge = 1
                    for vertex in tetEdges[edge].vertices:
                        tetVertices[vertex].isHullVertex = 1
                break
        
        if not tetHasHullFace:
            print("never found hull face")

fullConvexhull = ConvexHull(tetVerticesPositions)

hullFaces = []

for simp in fullConvexhull.simplices:
    hullFaces += [3, simp[0], simp[1], simp[2]] 
    

hullMesh = pv.PolyData()   
hullMesh.points = fullConvexhull.points
hullMesh.faces = np.array(hullFaces)   

#pl.add_mesh(hullMesh, color='manganeseblue')



points_poly = pv.PolyData(np.array(voronoiVertices))
select = points_poly.select_enclosed_points(hullMesh, check_surface=True, inside_out=True, tolerance=0.001)

pts = points_poly.extract_points(select['SelectedPoints'].view(bool), adjacent_cells=False)
#pl.add_points(pts, color='r')

for i in range(0, len(tetrahedra)):
    if select["SelectedPoints"][i]:

        crcm = tetrahedra[i].circumcenter #pir, index, position

        neigh = mesh.cell_neighbors(i, "faces") 

        accumulated = np.array([0.0, 0.0, 0.0])
        c = 0
        for n in neigh:
            ntet = tetrahedra[n]
            #if not ntet.isHullTet: continue
            c+=1
            accumulated += ntet.circumcenter[1]
        
        newCrcm = accumulated / c

        voronoiVertices[i] = newCrcm


# prints only surface faces
for face in tetFaces:
    if not face.isHullFace: continue
    tempFace = pv.PolyData()
    tempFace.points = face.points
    tempFace.faces = [3, 0, 1, 2]

    #pl.add_mesh(tempFace)



faceObjects = []
usedEdges = [0] * len(tetEdges)

print("pass making faces")

for i in range(0, len(tetrahedra)):

    tet = tetrahedra[i]

    
    for faceIndex in tet.faces:
        face = tetFaces[faceIndex]              
        for edgeIndex in face.edges:            
            if usedEdges[edgeIndex]: continue   #already used in face

            edge = tetEdges[edgeIndex]
            edgeTets = edge.tets
            vertices = []
            for edgeTet in edgeTets:
                # if isinstance(tetrahedra[edgeTet].hullFace, tuple):
                #     print("aqui")
                vertices.append(tetrahedra[edgeTet].circumcenter)    #esto es un par, con el indice y la posicion

            faceObjects.append((len(edgeTets), vertices)) 
            
            usedEdges[edgeIndex] = 1




testMesh = pv.PolyData()

voronoiFaces = []

print("going trough faces ordening tets so they follow convex hull topology")
for face in faceObjects:
    numOfVertices = face[0]
    faceVertices = face[1]

    if numOfVertices >= 4:
        positions = []

        for _, coords in faceVertices: positions.append(coords)

        hull = ConvexHull(positions, qhull_options="QJ")

        indices = []
        for simp in hull.simplices:
            simp = list(map(lambda x: faceVertices[x][0], simp))
            indices += [3] + simp
        

        voronoiFaces += indices

    if numOfVertices == 3:
        voronoiFaces += [3]
        for vertexIndex, coords in faceVertices: voronoiFaces += [vertexIndex]

    # # aqui tengo 2 circumcentros y quiero generar una cara con esos + los centros de cada una de las 2 cara superficial
    if numOfVertices == 2:
        positions = []


        in0 = []

        for i in range(2):
            positions.append(faceVertices[i][1]) # circumcentros ya conocidos
            # obtener centros de caras superficiales
            t0 = tetrahedra[faceVertices[i][0]]
            centerOfHullFace = tetFaces[t0.hullFace].center
            positions.append(centerOfHullFace)
            faceVertices.append([len(voronoiVertices), centerOfHullFace])
            voronoiVertices.append(centerOfHullFace)


        hull = ConvexHull(positions, qhull_options="QJ")

        indices = []
        for simp in hull.simplices:
            simp = list(map(lambda x: faceVertices[x][0], simp))
            indices += [3] + simp
        

        voronoiFaces += indices

    # if numOfVertices == 1:
    #     print("1!")
        

def findVoronoiVert(point, pointList):
    res = -1
    for i in range(0, len(pointList)):
        vpos = pointList[i]
        if abs(point[0] - vpos[0]) < 1e-9 and abs(point[1] - vpos[1]) < 1e-9 and abs(point[2] - vpos[2]) < 1e-9 : 
            res = i 
            break 
    return res


# quiero ir face por face, solo por las hull faces
# luego ir edge por edge, si no he visto ese edge obtengo la otra cara y hago la cara con los circumcentros de cada tet de las caras
seenEdge = [0] * len(tetEdges)
for i in range(0, len(tetFaces)):
    face = tetFaces[i]

    if face.isHullFace == 0: continue
    faceTet = tetrahedra[face.tets[0]]

    if not tetrahedra[face.tets[0]].isHullTet: faceTet = tetrahedra[face.tets[1]]

    for edge in face.edges:
        if seenEdge[edge]: continue

        edgeFaces = tetEdges[edge].faces
        pairFaceIndex = -1
        for edgeFace in edgeFaces:
            if edgeFace == i: continue 
            if tetFaces[edgeFace].isHullFace: 
                pairFaceIndex = edgeFace
                break 
            pairFaceIndex = edgeFace

    

        pairFace = tetFaces[pairFaceIndex]
        pairTet = tetrahedra[pairFace.tets[0]]
        if not tetrahedra[pairFace.tets[0]].isHullTet: pairTet = tetrahedra[pairFace.tets[1]]

        # add new vertices to voronoi vertices array
        newVert0 = findVoronoiVert(face.center, voronoiVertices)
        if newVert0 == -1: #nuevo vertice
            newVert0 = len(voronoiVertices)
            voronoiVertices.append(face.center)
        
        newVert1 = findVoronoiVert(pairFace.center, voronoiVertices)
        if newVert1 == -1:
            newVert1 = len(voronoiVertices)
            voronoiVertices.append(pairFace.center)


        pointIndices = [faceTet.circumcenter[0], pairTet.circumcenter[0], newVert0, newVert1]
        pointPositions = [faceTet.circumcenter[1], pairTet.circumcenter[1], face.center, pairFace.center]

        hull = ConvexHull(pointPositions, qhull_options="QJ")
        
        indices = []
        for simp in hull.simplices:
            simp = list(map(lambda x: pointIndices[x], simp))
            indices += [3] + simp
        

        voronoiFaces += indices

        seenEdge[edge] = 1


#suface faces
for vert in tetVertices:
    if not vert.isHullVertex: continue

    posiciones = []
    indices = []
    for vertFace in vert.faces:
        if not tetFaces[vertFace].isHullFace: continue
        posiciones.append(tetFaces[vertFace].center)
        indices.append(findVoronoiVert(tetFaces[vertFace].center, voronoiVertices))

    hullResult = []

    hull = ConvexHull(posiciones, qhull_options="QJ")
        
    for simp in hull.simplices:
        simp = list(map(lambda x: indices[x], simp))
        hullResult += [3] + simp
    

    #voronoiFaces += hullResult





        

# fullConvexhull = ConvexHull(tetVerticesPositions)

# hullFaces = []

# for simp in fullConvexhull.simplices:
#     hullFaces += [3, simp[0], simp[1], simp[2]] 
    

# hullMesh = pv.PolyData()   
# hullMesh.points = fullConvexhull.points
# hullMesh.faces = np.array(hullFaces)   



polyMesh = pv.PolyData()
polyMeshPoints = []
faceIndices = []
for face in tetFaces:

    vertices = face.points
    faceIndices += [3, len(polyMeshPoints), len(polyMeshPoints) + 1, len(polyMeshPoints) + 2]
    for vert in vertices:
        polyMeshPoints.append(vert)
    
polyMesh.points = np.array(polyMeshPoints)
polyMesh.faces = np.array(faceIndices)

#pl.add_mesh(polyMesh, color = 'manganeseblue')


print("finished")



voronoiMesh = pv.PolyData()
voronoiMesh.points = np.array(voronoiVertices)

voronoiMesh.faces = np.array(voronoiFaces) 


#pl.add_mesh(pts,color='r')

pl.add_mesh(voronoiMesh)
#pl.add_mesh(mesh, style="wireframe")
#pl.add_points(np.array(newborderVertices), render_points_as_spheres=True)
#pl.add_mesh(testMesh )

# pl.add_mesh(surfaceMesh,
#             show_edges=True)

pl.show()
