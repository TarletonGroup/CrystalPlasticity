# Nicolo Grilli
# University of Oxford
# AWE project 2020
# 4 giugno 2020

# main function

import numpy as np
from numpy import genfromtxt
from abaqusparser import AbaqusParser
from mesh import Mesh
from interfnode import Interfnode
from newmesh import NewMesh

# parse abaqus input file
inpfile = AbaqusParser('Job-1')
inpfile.ReadElsets()
inpfile.ReadNodes()
inpfile.ReadBulkElems()

# create mesh and parse file
m = Mesh('Job-1')
m.ReadElsets('Job-1')

# operations on the mesh
m.CreateConnectivity()
m.FindInterfNodes()
m.FindInterfElems()
m.CreateFaces()
m.CreateOrderCohNodes()

# create new mesh and operations
nm = NewMesh(m)
nm.CreateTwinNodes()
nm.MakeNewBulkElems()
nm.MakeNewInterfElems()

# write new mesh to file
nm.WriteNewNodeFile('Job-1')
nm.WriteNewElemsFile('Job-1')

# cohesive elements
nm.CreateCohesiveElems()
nm.WriteCohElems('Job-1')

# dummy elements
nm.WriteDummyElems('Job-1')









