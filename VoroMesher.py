import sys
import numpy as np
import pyvista as pv
from scipy.spatial import ConvexHull, Delaunay
from tqdm import tqdm

## UTILS FUNCTIONS

# Calculates a tetrahedron circumcenter from an array of vertices
# Implementation based on http://mathworld.wolfram.com/Circumsphere.html
def calc_tetrahedron_cirucmcenter(tet):

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

# Searches a list of Vertex for a point
def findPoint(point, verticesList):
    res = -1
    for i in range(0, len(verticesList)):
        vpos = verticesList[i].position
        if abs(point[0] - vpos[0]) < 1e-9 and abs(point[1] - vpos[1]) < 1e-9 and abs(point[2] - vpos[2]) < 1e-9 : 
            res = i 
            break 
    return res

# Searches a list of Edges for an adge that has the same center as edgeObject
def findEdge(edgeObject, edgesList):
    res = -1
    center = edgeObject.center
    for i in range(0, len(edgesList)):
        ecenter = edgesList[i].center
        if abs(center[0] - ecenter[0]) < 1e-9 and abs(center[1] - ecenter[1]) < 1e-9 and abs(center[2] - ecenter[2]) < 1e-9: 
            res = i 
            break 
    return res

# Searches a list of Faces for a face that has the same center as faceObject
def findFace(faceObject, facesList):
    res = -1
    center = faceObject.center
    for i in range(0, len(facesList)):
        fcenter = facesList[i].center
        if abs(center[0] - fcenter[0]) < 1e-9 and abs(center[1] - fcenter[1]) < 1e-9 and abs(center[2] - fcenter[2]) < 1e-9:
            res = i 
            break 
    return res




## MESH STRUCTURE CLASSES

# Represents a Delaunay Vertex
class Vertex:

    def __init__(self, p):
        self.position = p
        self.tets = []
        self.faces = []
        self.edges = []
        self.isHullVertex = 0

# Represents a Delaunay Edge
class Edge:

    def __init__(self, edgeObject):
        self.points = edgeObject.points
        self.center = edgeObject.center
        self.vertices = []
        self.faces = []
        self.tets = []
        self.isHullEdge = 0

# Represents a Delaunay Face
class Face:

    def __init__(self, faceObject):
        self.points = faceObject.points
        self.center = faceObject.center

        self.edges = []
        self.tets = []

        self.isHullFace = 0

# Represents a Tetrahedron
class Tetrahedro:



    def __init__(self, p, c):
        self.points = p
        self.center = c
        self.faces = []

        self.edges = []
        self.circumcenter = -1
        self.isHullTet = 0
        self.hullFace = -1

    def getCircumcenter(self):
        return calc_tetrahedron_cirucmcenter(self.points)



# Prints hull faces from list of faces
def renderHullFaces(plotter, facesList):
    for face in facesList:
        if not face.isHullFace: continue
        tempFace = pv.PolyData()
        tempFace.points = face.points
        tempFace.faces = [3, 0, 1, 2]

        plotter.add_mesh(tempFace)


def findVoronoiVert(point, pointList):
    res = -1
    for i in range(0, len(pointList)):
        vpos = pointList[i]
        if abs(point[0] - vpos[0]) < 1e-9 and abs(point[1] - vpos[1]) < 1e-9 and abs(point[2] - vpos[2]) < 1e-9 : 
            res = i 
            break 
    return res


def initializeStructureArrays(pyvistaMesh, listOfTets, listOfFaces, listOfEdges, listOfVertices, listOfVerticesPositions, listOfVoronoiVertices):
    #print("VoroMesher: Initializing structure arrays")

    for i in tqdm(range(pyvistaMesh.n_cells), desc="VoroMesher: Initializing structure arrays", leave=False):

    #for i in range(0, pyvistaMesh.n_cells):

        # get cell
        cell = pyvistaMesh.get_cell(i)

        # make new tetrahedra
        tet = Tetrahedro(cell.points, cell.center)

        # iterate faces
        cell_faces = pyvistaMesh.get_cell(i).faces
        for face in cell_faces:

            # make new face
            f0 = Face(face)

            newFace = 0
            # search if face exists already
            index = findFace(f0, listOfFaces)
            if index == -1: # doesnt exists, gets new index
                newFace = 1
                index = len(listOfFaces) #apendiar al final
            else:           # exists, gets face from array
                f0 = listOfFaces[index]


            edgeObjects = []
            # get face edges
            faceEdges = face.edges
            # iterate face edges 
            for edge in faceEdges:
                # new edge
                e0 = Edge(edge)

                # check if edge exists
                newEdge = 0
                edgeIndex = findEdge(e0, listOfEdges)
                if edgeIndex == -1:
                    newEdge = 1
                    edgeIndex = len(listOfEdges)
                else:
                    e0 = listOfEdges[edgeIndex]


                vertexObjects = []
                # get edge vertices
                edgeVertices = edge.points
                # iterate vertices
                for vertex in edgeVertices:
                    # make new vertex
                    p0 = Vertex(vertex) #vertex is array of 3 points

                    # check if it exists
                    point0Index = findPoint(vertex, listOfVertices)
                    if point0Index == -1:
                        point0Index = len(listOfVertices)
                        p0.tets.append(i)
                        p0.faces.append(index)
                        p0.edges.append(edgeIndex)
                        listOfVertices.append(p0)
                        listOfVerticesPositions.append(np.array(vertex))
                    else:
                        p0 = listOfVertices[point0Index]
                        if not (i in p0.tets): p0.tets.append(i)
                        if not(index in p0.faces): p0.faces.append(index)
                        if not(edgeIndex in p0.edges): p0.edges.append(edgeIndex)

                    vertexObjects.append(point0Index)
                # add vertices to edge
                for vertex in vertexObjects:
                    if not (vertex in e0.vertices): e0.vertices.append(vertex) 

                # add face to edge
                if not (index in e0.faces) : e0.faces.append(index)     
                # add tet to edge
                if not (i in e0.tets) : e0.tets.append(i)               

                # add edge to edge list
                if newEdge: listOfEdges.append(e0)

                edgeObjects.append(edgeIndex)

            # add edges to face
            for edge in edgeObjects:
                if not (edge in f0.edges): f0.edges.append(edge) #add edges to face

            # add tet to face
            if not (i in f0.tets): f0.tets.append(i)        # add tet to face

            # add face to faces list
            if newFace: listOfFaces.append(f0)             # se agrega a la lista de faces si es que no existe

            # add face to tet
            tet.faces.append(index)


        # gets tet circumcenter
        tetCircum = tet.getCircumcenter()
        # save circumcenter with index and position on tet
        tet.circumcenter = (len(listOfVoronoiVertices), tetCircum)
        # save position to list
        listOfVoronoiVertices.append(tetCircum)

        # add tetrahedra
        listOfTets.append(tet)


    #print("VoroMesher: Marking surface tetrahedra and faces")
    # Identifying and marking surface tets and faces
    for i in tqdm(range(pyvistaMesh.n_cells), desc="VoroMesher: Marking surface tetrahedra and faces", leave=False):
    #for i in range(0, pyvistaMesh.n_cells):                            # for each tertahedra
        neigh = pyvistaMesh.cell_neighbors(i, "faces")                 # get neighbors
        if len(neigh) < 4:                                              # if a tetrahedra has less than 4 neighbors, its a border tet

            # marking tet
            listOfTets[i].isHullTet = 1
    
            # Iterate faces, must find hull face out of the 4 faces
            for face in listOfTets[i].faces:     
                # mark as hull face, until is not
                found = 0                                      

                # for each neighbor
                for neighbor in neigh: 
                    # get neighbor faces
                    neighFaces = listOfTets[neighbor].faces     
                    # if i found a face in the neighbor faces, its not the hull face
                    if face in neighFaces:  
                        found = 1
                        break
                # if i went trough all neighbors without finding the face, its the hull face
                if not found:
                    listOfFaces[face].isHullFace = 1

                    # mark edges and vertices as hull too
                    for edge in listOfFaces[face].edges: 
                        listOfEdges[edge].isHullEdge = 1
                        for vertex in listOfEdges[edge].vertices:
                            listOfVertices[vertex].isHullVertex = 1

                    if listOfTets[i].hullFace != -1:
                        tetPair = (listOfTets[i].hullFace, face)
                        listOfTets[i].hullFace = tetPair
                        break
                    else:
                        listOfTets[i].hullFace = face

def suppressWildCircumcenters(listOfTets, listOfVerticesPositions, listOfVoronoiVertices):
    #print("VoroMesher: Suppressing circumcenters outside of Convex Hull")

    # Calculates the convex hull of original mesh, 
    # we use this to identify circumcenters too far from the mesh (outisde of vertex hull)
    fullConvexhull = ConvexHull(listOfVerticesPositions)

    # holds surface faces on this format [x, i, j, k] where x is the amount of vertices (always 3 in this case), and i, j, k the vertex's indices
    hullFaces = []

    # defines faces with convex hull
    for simp in fullConvexhull.simplices:
        hullFaces += [3, simp[0], simp[1], simp[2]] 
        
    # makes a mesh to show surface
    hullMesh = pv.PolyData()   
    hullMesh.points = fullConvexhull.points
    hullMesh.faces = np.array(hullFaces) 

    # Gets and shows voronoi vertices outside of Convex Hull
    points_poly = pv.PolyData(np.array(listOfVoronoiVertices))
    select = points_poly.select_enclosed_points(hullMesh, check_surface=True, inside_out=True, tolerance=0.001)

    for i in tqdm(range(len(listOfTets)), desc="VoroMesher: Suppressing circumcenters outside of Convex Hull", leave=False):
    #for i in range(0, len(listOfTets)):
        if select["SelectedPoints"][i]:
            # sets position as average position from neighboring voronoi vertices
            listOfVoronoiVertices[i] = listOfTets[i].center

    pts = points_poly.extract_points(select['SelectedPoints'].view(bool), adjacent_cells=False)

    return pts

def makeInternalVoronoiFaces(listOfTets, listOfFaces, listOfEdges, listOfVoronoiVertices, listOfVoronoiFaces):
    #print("VoroMesher: Defining internal Voronoi face vertices")

    # list of tuples(n of vertices, list(vertices))
    faceObjects = []
    # list to set used edges
    usedEdges = [0] * len(listOfEdges)

    for i in tqdm(range(len(listOfTets)), desc="VoroMesher: Defining internal Voronoi face vertices", leave=False):
    #for i in range(0, len(listOfTets)): # iterate tetrahedra

        tet = listOfTets[i]

        # iterate tet faces
        for faceIndex in tet.faces:
            face = listOfFaces[faceIndex]   
            # iterate face edges           
            for edgeIndex in face.edges:            
                if usedEdges[edgeIndex]: continue   #there is only 1 face per edge

                edge = listOfEdges[edgeIndex]
                # get tets surrounding edge
                edgeTets = edge.tets
                vertices = []
                for edgeTet in edgeTets:
                    t0 = listOfTets[edgeTet]
                    vertices.append(t0.circumcenter)    # add tet circumcenter to the face's vertices
                    if edge.isHullEdge and t0.isHullTet:
                        if isinstance(t0.hullFace, tuple):
                            vertices.append((len(listOfVoronoiVertices), listOfFaces[t0.hullFace[0]].center))
                            listOfVoronoiVertices.append(listOfFaces[t0.hullFace[0]].center)

                            vertices.append((len(listOfVoronoiVertices), listOfFaces[t0.hullFace[1]].center))
                            listOfVoronoiVertices.append(listOfFaces[t0.hullFace[1]].center)
                        else:
                            vertices.append((len(listOfVoronoiVertices), listOfFaces[t0.hullFace].center))
                            listOfVoronoiVertices.append(listOfFaces[t0.hullFace].center)

                faceObjects.append((len(edgeTets), vertices)) # add face to list
                
                usedEdges[edgeIndex] = 1


    #print("VoroMesher: Making triangles for a Voronoi face")
    # iterate list of faces (set of vertices)
    for i in tqdm(range(len(faceObjects)), desc="VoroMesher: Making triangles for a Voronoi face", leave=False):
        face = faceObjects[i]
    #for face in faceObjects:
        numOfVertices = face[0]
        faceVertices = face[1]

        # regular case, face is in the inside of the mesh
        if numOfVertices >= 1:
            positions = []
            voronoiIndices = []

            nofhulltets = 0
            for vertexIndex, coords in faceVertices: 
                positions.append(coords)
            
            if len(faceVertices) == 3: 
                listOfVoronoiFaces += [3, faceVertices[0][0], faceVertices[1][0], faceVertices[2][0]]
                continue

            hull = ConvexHull(positions, qhull_options="QJ")

            indices = []
            for simp in hull.simplices:
                simp = list(map(lambda x: faceVertices[x][0], simp))
                indices += [3] + simp
            
            

            listOfVoronoiFaces += indices
            

#suface faces
def makeExternalVoronoiFaces(listOfVertices, listOfFaces, listOfVoronoiVertices, listOfVoronoiFaces):
    #print("VoroMesher: Making external Voronoi faces")

    for i in tqdm(range(len(listOfVertices)), desc="VoroMesher: Making external Voronoi faces", leave=False):
        vert = listOfVertices[i]
    #for vert in listOfVertices:
        if not vert.isHullVertex: continue

        posiciones = []
        indices = []
        for vertFace in vert.faces:
            if not listOfFaces[vertFace].isHullFace: continue
            posiciones.append(listOfFaces[vertFace].center)
            indices.append(findVoronoiVert(listOfFaces[vertFace].center, listOfVoronoiVertices))

        if len(posiciones) <= 3:
            if len(posiciones) == 3:
                listOfVoronoiFaces += [3, indices[0], indices[1], indices[2]]
            continue

        hullResult = []

        hull = ConvexHull(posiciones, qhull_options="QJ")
            
        for simp in hull.simplices:
            simp = list(map(lambda x: indices[x], simp))
            hullResult += [3] + simp
        

        listOfVoronoiFaces += hullResult

def main():

    if '-h' in sys.argv:
        print("Usage: [VTK file] [options ...]\n" \
        "where options include:\n" \
        "-h         print this message\n" \
        "-o         shows wireframe of original mesh\n" \
        "-ne        doesn't show external voronoi faces\n" \
        "-s         suppresses circumcenters outside of convex hull\n" \
        "-p         if -p, shows those circumcenters\n")

        sys.exit()

    pl = pv.Plotter()

    # Read a vtk file pass as argument to a mesh
    mesh = pv.read(sys.argv[1])

    if '-o' in sys.argv:
        pl.add_mesh(mesh, style="wireframe")

    # Initialize array of structures
    tetrahedra = []
    tetFaces = []
    tetEdges = []

    tetVertices = []
    tetVerticesPositions = []

    voronoiVertices = []
    voronoiFaces = []

    initializeStructureArrays(mesh, tetrahedra, tetFaces, tetEdges, tetVertices, tetVerticesPositions, voronoiVertices)

    if '-s' in sys.argv:
        pts = suppressWildCircumcenters(tetrahedra, tetVerticesPositions, voronoiVertices)      

        if '-p' in sys.argv:
            pl.add_points(pts, color='r')

    makeInternalVoronoiFaces(tetrahedra, tetFaces, tetEdges, voronoiVertices, voronoiFaces)

    if not '-ne' in sys.argv:
        makeExternalVoronoiFaces(tetVertices, tetFaces, voronoiVertices, voronoiFaces)


    voronoiMesh = pv.PolyData()
    voronoiMesh.points = np.array(voronoiVertices)
    voronoiMesh.faces = np.array(voronoiFaces) 

    pl.add_mesh(voronoiMesh)
    pl.show()



main()