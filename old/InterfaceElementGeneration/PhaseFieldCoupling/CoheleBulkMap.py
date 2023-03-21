# Nicolo Grilli
# University of Oxford
# AWE project 2020
# 28 Luglio 2020

# build a lookup table
# to find the neighbouring elements
# of a cohesive element
# as a fortran matrix to include in the UEL

import numpy as np
from numpy import genfromtxt

nodefilename = 'Square20umEl0p25-node.inp'
bulkelemfilename = 'Square20umEl0p25-bulk-elems.inp'
intelemfilename = 'Square20umEl0p25-int-elems-new.inp'

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

# the lookup table has:
# 1st column: index of interface element
# 2nd column: index of the first neighbour bulk element
# 3rd column: index of the second neighbour bulk element
cohelebulkmap = np.zeros(shape=(len(intelem),3))
for indiceintelem in range(0,len(intelem)):
	tempintnodes = intelem[indiceintelem,1:9]
	temp2neighbours = FindTwoNeighbourBulkElem(tempintnodes,bulkelem)
	cohelebulkmap[indiceintelem,0] = intelem[indiceintelem,0]
	cohelebulkmap[indiceintelem,1] = temp2neighbours[0]+1 # Abaqus index
	cohelebulkmap[indiceintelem,2] = temp2neighbours[1]+1 # Abaqus index

# write to file the lookup table
# in fortran matrix format
fid = open('CoheleBulkMap.f','w')

for indiceintelem in range(0,len(intelem)):
	fid.write('      ') # 6 empty columns for fortran 77
	fid.write('cohelebulkmap('+str(int(indiceintelem+1))+',:) = (/')
	fid.write(str(int(cohelebulkmap[indiceintelem,0])))
	fid.write(',')
	fid.write(str(int(cohelebulkmap[indiceintelem,1])))
	fid.write(',')
	fid.write(str(int(cohelebulkmap[indiceintelem,2])))
	fid.write('/)')
	fid.write('\n')

fid.close()











