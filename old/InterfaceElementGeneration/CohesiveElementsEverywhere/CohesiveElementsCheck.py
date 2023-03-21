# Nicolo Grilli
# University of Oxford
# 15 Maggio 2020

# check that orientation of the normal
# in the cohesive elements (GB) is correct

import numpy as np
from numpy import genfromtxt

nodefilename = 'Square20umEl0p25-node.inp'
bulkelemfilename = 'Square20umEl0p25-bulk-elems.inp'
intelemfilename = 'Square20umEl0p25-int-elems.inp'

nodes = genfromtxt(nodefilename,dtype=float,delimiter=",",skip_header=0,skip_footer=0)
bulkelem = genfromtxt(bulkelemfilename,dtype=int,delimiter=",",skip_header=0,skip_footer=0)
intelem = genfromtxt(intelemfilename,dtype=int,delimiter=",",skip_header=0,skip_footer=0)

# function to check if four nodes of a cohesive element
# correspond to four nodes of a bulk element
# bulknodes is a 8x1 numpy array indicating the nodes in one bulk element
# intnodes is a 8x1 numpy array indicating the nodes in one interface
def AreNodesInBulkElem(intnodes,bulknodes):
	count = 0 # count how many nodes are also in bulk elem
	for x in range(len(intnodes)):
		if intnodes[x] in bulknodes: # this node is also in bulk elem
			count = count + 1
	if (count == 4):
		return True
	else:
		return False

# function to
# find the two bulk elements that have four common nodes
# with a specific cohesive element
# intnodes is a 8x1 numpy array indicating the nodes in one interface
# bulkelem are all the bulk elements (indices and nodes)
def FindTwoNeighbourBulkElem(intnodes,bulkelem):
	TwoNeighbourBulkElem = np.zeros(shape=(2))
	count = 0 # index of TwoNeighbourBulkElem
	for indicebulkelem in range(0,len(bulkelem)):
		a = bulkelem[indicebulkelem,1:9]
		if (AreNodesInBulkElem(intnodes,a)):
			TwoNeighbourBulkElem[count] = indicebulkelem
			count = count + 1
	return TwoNeighbourBulkElem

# find the normal of a cohesive element
# assuming right hand rule
# moving in a circle along nodes 1,2,3,4 of
# that cohesive element
# intnodes is a 8x1 numpy array indicating the nodes in one interface
# nodes are all the nodes with coordinates
def NormalCohesiveElem(intnodes,nodes):
	# define local frame
	# xlocal = x axis from node 1 to node 2 of int elem
	# ylocal = y axis from node 1 to node 4 of int elem
	xlocal = np.zeros(shape=(3))
	ylocal = np.zeros(shape=(3))
	node1 = np.zeros(shape=(3))
	node2 = np.zeros(shape=(3))
	node4 = np.zeros(shape=(3))
	node1 = nodes[intnodes[0]-1,1:4] # coordinates of node 1
	node2 = nodes[intnodes[1]-1,1:4] # coordinates of node 2
	node4 = nodes[intnodes[3]-1,1:4] # coordinates of node 4
	xlocal = node2 - node1
	ylocal = node4 - node1
	zlocal = np.cross(xlocal,ylocal)	
	return zlocal

# check of the normal of the cohesive element
# points towards the second element
# which is the ones that share with the cohesive element
# the last four nodes of the cohesive element
# intnodes is a 8x1 numpy array indicating the nodes in one interface
# nodes are all the nodes with coordinates
# bulkelem are all the bulk elements (indices and nodes)
def IsNormalRightHand(intnodes,nodes,bulkelem):
	# find indices of the two neighbouring elements
	TwoNeighbourBulkElem = FindTwoNeighbourBulkElem(intnodes,bulkelem)
	# calculate normal of the interface element
	# assuming right hand rule
	# moving in a circle along nodes 1,2,3,4 of
	# that cohesive element
	zlocal = NormalCohesiveElem(intnodes,nodes)
	# which of the two elements share with the interface element
	# the first 4 nodes of the interface element?
	# the normal to the interface element should point on
	# the opposite direction with respect to that element
	indicefirstelem = int(TwoNeighbourBulkElem[0])
	indicesecondelem = int(TwoNeighbourBulkElem[1])
	tempfirstbulkelem = bulkelem[indicefirstelem,1:9]
	tempsecondbulkelem = bulkelem[indicesecondelem,1:9]
	firstnodeintelem = intnodes[0]
	if firstnodeintelem in tempfirstbulkelem: # check first neighbour
		tempbulkelem = tempfirstbulkelem
	if firstnodeintelem in tempsecondbulkelem: # check second neighbour
		tempbulkelem = tempsecondbulkelem
	# search node of the bulk element that does not belong to interface element
	foundnode = False
	indiceNotBelong = 0
	for x in range(0,8):	
		if (foundnode == False):
			if tempbulkelem[x] not in intnodes:
				indiceNotBelong = x
				foundnode = True
	fifthnodeintelem = intnodes[4]
	# construct vector from node 1 of interface element 
	# to node of the bulk element that does not belong to interface
	node1int = np.zeros(shape=(3))
	node5bulk = np.zeros(shape=(3))
	node1int = nodes[firstnodeintelem-1,1:4]
	node5bulk = nodes[tempbulkelem[indiceNotBelong]-1,1:4]
	zinttobulk = node5bulk - node1int
	# zinttobulk should point in the opposite direction
	# with respect to the interface element normal
	if (zlocal.dot(zinttobulk) > 0):
		return False
	else:
		return True

# function to change the orientation of an interface element
# surface from left hand rule to right hand rule
# intnodes is a 8x1 numpy array indicating the nodes in one interface
def ChangeIntElemOrientation(intnodes):
	newintnodes = np.zeros(shape=(8))
	for x in range(0,2):
		xper4 = 4*x
		newintnodes[0+xper4] = intnodes[0+xper4]
		newintnodes[1+xper4] = intnodes[3+xper4]
		newintnodes[2+xper4] = intnodes[2+xper4]
		newintnodes[3+xper4] = intnodes[1+xper4]
	return newintnodes

# check for which interface elements the normal
# points in the wrong direction and change it
newintelem = np.zeros(shape=(len(intelem),9))
for indiceintelem in range(0,len(intelem)):
	tempintnodes = intelem[indiceintelem,1:9]
	if (IsNormalRightHand(tempintnodes,nodes,bulkelem) == False):
		print(indiceintelem)
		newintnodes = ChangeIntElemOrientation(tempintnodes)
		newintelem[indiceintelem,1:9] = newintnodes[0:8]
		newintelem[indiceintelem,0] = intelem[indiceintelem,0]
	else:
		newintelem[indiceintelem,0:9] = intelem[indiceintelem,0:9]

# write to file new interface elements
np.savetxt('Job-1-int-elem-new.inp', newintelem, fmt='%d', delimiter=',')







			 


