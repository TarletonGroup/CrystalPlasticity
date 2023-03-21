# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2017 replay file
# Internal Version: 2016_09_27-22.54.59 126836
# Run by engs1992 on Tue Jun 16 14:20:49 2020
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
import numpy as np
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=218.485397338867, 
    height=100.585189819336)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
Mdb()

# parse geometrical parameters from file
fparam = open("Parameters.txt","r")

pwidth = float(fparam.readline())
pheight = float(fparam.readline())
pdepth = float(fparam.readline())
Ngrains = int(fparam.readline())
caepath = fparam.readline()

fparam.close()

#: A new model database has been created.
#: The model "Model-1" has been created.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)

# draw rectangle
s.rectangle(point1=(0.0, 0.0), point2=(pwidth,pheight))

p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-1']

# extrude to create parallelepiped
p.BaseSolidExtrude(sketch=s, depth=pdepth)

s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']

p = mdb.models['Model-1'].parts['Part-1']

# name of files with Euler angles, vertices and edges
euleranglesfilename = "n" + str(Ngrains) + "-id1-eulerangles.tess"
datumpointsfilename = "n" + str(Ngrains) + "-id1-datumpoints.tess"
edgesfilename = "n" + str(Ngrains) + "-id1-edges.tess"

# read external datum points file and add datum points to cae model
fid = open(datumpointsfilename,"r")

# arrays of x, y, z coordinates of vertices
arrx = np.zeros(shape=(0)) 
arry = np.zeros(shape=(0))

count = 1
for line in fid:
	tempdata = line.split()
	# get x and y coordinate of datum point
	tempx = float(tempdata[1])
	tempy = float(tempdata[2])
	arrx = np.append(arrx,tempx)
	arry = np.append(arry,tempy)
	p.DatumPointByCoordinate(coords=(tempx, tempy, pdepth))
	count = count + 1

fid.close()

p = mdb.models['Model-1'].parts['Part-1']
f = p.faces

# define boundingbox to pick the entire upper face
BoundBox1 = np.zeros(shape=(3))
BoundBox2 = np.zeros(shape=(3))

BoundBox1[0] = -0.1*pwidth
BoundBox1[1] = -0.1*pheight
BoundBox1[2] = 0.9*pdepth

BoundBox2[0] = 1.1*pwidth
BoundBox2[1] = 1.1*pheight
BoundBox2[2] = 1.1*pdepth

pickedFaces = f.getByBoundingBox(BoundBox1[0],BoundBox1[1],BoundBox1[2],BoundBox2[0],BoundBox2[1],BoundBox2[2])
v, e, d = p.vertices, p.edges, p.datums

# array that contains point at the boundary of edges
# using .tess file numbering
edgepoint1 = np.array([],dtype='uint64')
edgepoint2 = np.array([],dtype='uint64')
# array with the edge indices
# counting only edges defined on the upper surface
# and leaving 0 the ones on the corner edges
edgeindex = np.array([],dtype='uint64')

fied = open(edgesfilename,"r")

count = 1

for line in fied:
	InnerEdge = True # flag to indicate that this is not an edge at the boundary of the surface 
	p = mdb.models['Model-1'].parts['Part-1']
	f = p.faces
	pickedFaces = f.getByBoundingBox(BoundBox1[0],BoundBox1[1],BoundBox1[2],BoundBox2[0],BoundBox2[1],BoundBox2[2])
	v, e, d = p.vertices, p.edges, p.datums
	tempdata = line.split()
	# get indices of the datum points bounding the edge
	# indices start from 1
	# these are indices in the .tess file, not abaqus
	temppoint1 = int(tempdata[1]) # corresponds to abaqus datum point number - 1
	temppoint2 = int(tempdata[2]) # corresponds to abaqus datum point number - 1
	edgepoint1 = np.append(edgepoint1,np.uint64(temppoint1))
	edgepoint2 = np.append(edgepoint2,np.uint64(temppoint2))
	# is this an edge at the boundary?
	if (arrx[temppoint1-1] == arrx[temppoint2-1]):
		if (arrx[temppoint1-1] == 0.0 or arrx[temppoint1-1] == pwidth):
			InnerEdge = False
	if (arry[temppoint1-1] == arry[temppoint2-1]):
		if (arry[temppoint1-1] == 0.0 or arry[temppoint1-1] == pheight):
			InnerEdge = False
	if (InnerEdge):
		# find abaqus indices
		abqindex1 = temppoint1+1
		abqindex2 = temppoint2+1	
		p.PartitionFaceByShortestPath(point1=d[abqindex1], point2=d[abqindex2], 
			faces=pickedFaces)
		edgeindex = np.append(edgeindex,np.uint64(count))
		count = count + 1
	else:
		edgeindex = np.append(edgeindex,np.uint64(0))

fied.close()

# define boundingbox to pick the entire volume
BoundBox1[0] = -0.1*pwidth
BoundBox1[1] = -0.1*pheight
BoundBox1[2] = -0.1*pdepth

BoundBox2[0] = 1.1*pwidth
BoundBox2[1] = 1.1*pheight
BoundBox2[2] = 1.1*pdepth

p = mdb.models['Model-1'].parts['Part-1']
c = p.cells
pickedCells = c.getByBoundingBox(BoundBox1[0],BoundBox1[1],BoundBox1[2],BoundBox2[0],BoundBox2[1],BoundBox2[2])
alledges, d = p.edges, p.datums

# create grains = cells
count = 0
for indexe in range(0,len(edgeindex)): # cycle over all edges in tess file
	if (edgeindex[indexe] > 0): # corresponds to an edge on the upper surface
		temppoint1 = edgepoint1[indexe]
		temppoint2 = edgepoint2[indexe]
		midx = (arrx[temppoint1-1] + arrx[temppoint2-1])/2.0
		midy = (arry[temppoint1-1] + arry[temppoint2-1])/2.0
		p = mdb.models['Model-1'].parts['Part-1']
		c = p.cells
		pickedCells = c.getByBoundingBox(BoundBox1[0],BoundBox1[1],BoundBox1[2],BoundBox2[0],BoundBox2[1],BoundBox2[2])
		alledges, d = p.edges, p.datums
		# find current edge on the upper surface
		tempedge = alledges.findAt((midx,midy,pdepth))
		sweepEdge = alledges.findAt((pwidth,pheight,0.4*pdepth))
		p.PartitionCellByExtrudeEdge(line=sweepEdge, cells=pickedCells, 
			edges=tempedge, sense=REVERSE)

# read Euler angles
phi1 = np.zeros(shape=(0))
Phi = np.zeros(shape=(0))
phi2 = np.zeros(shape=(0))

feuler = open(euleranglesfilename,"r")

for line in feuler:
	tempdata = line.split()
	# get coordinates of the seed
	tempphi1 = float(tempdata[0])
	tempPhi = float(tempdata[1])
	tempphi2 = float(tempdata[2])
	phi1 = np.append(phi1,tempphi1)
	Phi = np.append(Phi,tempPhi)
	phi2 = np.append(phi2,tempphi2)

feuler.close()

# generate rotation matrix
# from Euler angles
# angles in radians
def bungeMatrix(phi1,Phi,phi2):
	matrix = np.zeros(shape=(3,3))
	matrix[0,0] = np.cos(phi1) * np.cos(phi2) - np.sin(phi1) * np.sin(phi2) * np.cos(Phi)
	matrix[0,1] = - np.cos(phi1) * np.sin(phi2) - np.sin(phi1) * np.cos(phi2) * np.cos(Phi)
	matrix[0,2] = np.sin(phi1) * np.sin(Phi)
	matrix[1,0] = np.sin(phi1) * np.cos(phi2) + np.cos(phi1) * np.sin(phi2) * np.cos(Phi)
	matrix[1,1] = - np.sin(phi1) * np.sin(phi2) + np.cos(phi1) * np.cos(phi2) * np.cos(Phi)
	matrix[1,2] = - np.cos(phi1) * np.sin(Phi)
	matrix[2,0] = np.sin(phi2) * np.sin(Phi)
	matrix[2,1] = np.cos(phi2) * np.sin(Phi)
	matrix[2,2] = np.cos(Phi)
	return matrix
			
# create materials and sections
# and assign them to cells

session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
    engineeringFeatures=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
	
p = mdb.models['Model-1'].parts['Part-1']
c = p.cells

tempgrainnumber = 0

for cella in p.cells:
	tempgrainnumber = tempgrainnumber + 1
	tempgrainname = 'GRAIN' + str(tempgrainnumber)
	mdb.models['Model-1'].Material(name=tempgrainname)
	mdb.models['Model-1'].materials[tempgrainname].Depvar(n=125)
	# assign rotation matrix
	tempphi1 = (np.pi/180.0) * phi1[tempgrainnumber-1]
	tempPhi = (np.pi/180.0) * Phi[tempgrainnumber-1]
	tempphi2 = (np.pi/180.0) * phi2[tempgrainnumber-1]
	temprot = bungeMatrix(tempphi1,tempPhi,tempphi2)
	mechConstTuple = (5,temprot[0,0],temprot[0,1],temprot[0,2],temprot[1,0],temprot[1,1],temprot[1,2], \
		temprot[2,0],temprot[2,1],temprot[2,2],tempgrainnumber)
	mdb.models['Model-1'].materials[tempgrainname].UserMaterial(mechanicalConstants=mechConstTuple)
	# create section
	mdb.models['Model-1'].HomogeneousSolidSection(name=tempgrainname, material=tempgrainname, 
	    thickness=None)
	p = mdb.models['Model-1'].parts['Part-1']
	c = p.cells
	puntocella = c[tempgrainnumber-1].pointOn
	cells = c.findAt(((puntocella[0][0],puntocella[0][1],puntocella[0][2]),))
	tempregion = p.Set(cells=cells, name=tempgrainname)
	region = regionToolset.Region(cells=cells)
	p = mdb.models['Model-1'].parts['Part-1']
	p.SectionAssignment(region=tempregion, sectionName=tempgrainname, offset=0.0, 
	    offsetType=MIDDLE_SURFACE, offsetField='', 
	    thicknessAssignment=FROM_SECTION)


