# Nicolo Grilli
# University of Oxford
# AWE project 2020
# 4 giugno 2020

# class representing the mesh

import numpy as np
from numpy import genfromtxt
from element import Element
from node import Node
from interfnode import Interfnode
from interfelement import InterfElement

class Mesh:

	# prefissofile is the prefix of nodes
	# and bulk elements files
	def __init__(self,prefissofile):
		# read files with nodes and elements
		nodefilename = prefissofile + '-node.inp'
		bulkelemfilename = prefissofile + '-bulk-elems.inp'
		nodefile = genfromtxt(nodefilename,dtype=float,delimiter=",",skip_header=0,skip_footer=0)
		bulkelemfile = genfromtxt(bulkelemfilename,dtype=int,delimiter=",",skip_header=0,skip_footer=0)
		# define a list of nodes and elements
		objnodes = []
		objbulkelem = []
		for linea in range(0,len(nodefile)):
			objnodes.append(Node(nodefile[linea,0],nodefile[linea,1:4]))
		for linea in range(0,len(bulkelemfile)):
			objbulkelem.append(Element(bulkelemfile[linea,0],bulkelemfile[linea,1:9]))
		self.nodes = objnodes
		self.bulkelements = objbulkelem
		self.interfnodes = []
		self.interfelems = []

	# assign grains to elements
	def ReadElsets(self,prefissofile):
		nomefile = prefissofile + '-elsets.inp'
		fid = open(nomefile)
		indicegrain = np.zeros(shape=(0))
		intervstart = np.zeros(shape=(0))
		intervend = np.zeros(shape=(0))
		for lineagrain in fid:
			if (lineagrain[0] == '*'):
				listakeywords = lineagrain.split(',')
				grainkeyword = listakeywords[1]
				posgrain = grainkeyword.find('GRAIN')
				indicegrain = np.append(indicegrain,int(grainkeyword[posgrain+5:]))
			else:
				eleminterval = lineagrain.split(',')
				intervstart = np.append(intervstart,int(eleminterval[0]))
				intervend = np.append(intervend,int(eleminterval[1]))
		fid.close()
		# check all elements and assign grain number 
		for bulkelem in self.bulkelements:
			for grain in range(0,len(indicegrain)):
				if (bulkelem.index >= intervstart[grain] and bulkelem.index <= intervend[grain]):
					bulkelem.grain = int(indicegrain[grain])
		return 1

	# assign to each node the elements
	# to which the node belongs
	def CreateConnectivity(self):
		for bulkelem in self.bulkelements: # cycle over bulk elements
			for indicenode in range(0,8):
				tempnode = bulkelem.nodes[indicenode] # tempnode is abaqus node index
				tempobjnode = self.nodes[tempnode-1]
				tempobjnode.connectivity = np.append(tempobjnode.connectivity,int(bulkelem.index))
		return 1

	# search through all nodes and find the ones
	# at the grain boundaries
	def FindInterfNodes(self):
		objinterfnodes = [] # list of interface nodes
		for node in self.nodes:
			tempindex = node.index
			tempcoords = node.coords
			tempconnectivity = np.int_(node.connectivity)
			nodegrains = np.zeros(shape=(len(tempconnectivity))) # grains of the neighbouring elements
			count = 0
			for elem in tempconnectivity: # elem is abaqus element index
				nodegrains[count] = int(self.bulkelements[elem-1].grain)
				count = count + 1
			# check if all grains in the elements surrounding the node are the same
			IsBulkNode = np.all(nodegrains == nodegrains[0])
			if (not IsBulkNode): # this is an interface node
				objinterfnodes.append(Interfnode(tempindex,tempcoords))
				objinterfnodes[-1].connectivity = tempconnectivity
				# calculate multiplicity, 3 if it is at a triple junction
				objinterfnodes[-1].multiplicity = len(np.unique(nodegrains))
				# this node will belong to the grain with smallest index
				objinterfnodes[-1].grain = int(min(nodegrains))
				# assign new connectivity
				count = 0
				for elem in tempconnectivity:
					if (nodegrains[count] == objinterfnodes[-1].grain):
						objinterfnodes[-1].newconnectivity = np.append(objinterfnodes[-1].newconnectivity,int(elem))
					count = count + 1
		self.interfnodes = objinterfnodes
		return 1

	# search through the elements and find the ones
	# that have four nodes on the grain boundaries
	# in common with another element
	# note that one element can have four nodes on a first GB
	# and four other nodes on a second GB
	# in this case it will be counted twice
	def FindInterfElems(self):
		objinterfelems = [] # list of interface elements
		indexobjinterfelems = [] # keep a list of indices of elements pairs identified to avoid repetition
		for node in self.interfnodes:
			tempconnectivity = np.int_(node.connectivity)
			# check all the pairs and see if two elements in different grains
			# share four nodes
			for x in tempconnectivity:
				xnodegrains = int(self.bulkelements[x-1].grain)
				for y in tempconnectivity:
					ynodegrains = int(self.bulkelements[y-1].grain)
					if (xnodegrains != ynodegrains): # excludes the same element case
						xnodes = self.bulkelements[x-1].nodes
						ynodes = self.bulkelements[y-1].nodes
						countcommonnodes = 0
						for c in range(0,8):
							if xnodes[c] in ynodes:
								countcommonnodes = countcommonnodes + 1
						if (countcommonnodes == 4): # this is a common face on the GB
							if ([x,y] not in indexobjinterfelems): # not already included
								if ([y,x] not in indexobjinterfelems): # not already included
									objinterfelems.append(InterfElement(x,xnodes))
									objinterfelems[-1].shareface = y
									objinterfelems[-1].grain = xnodegrains
									objinterfelems.append(InterfElement(y,ynodes))
									objinterfelems[-1].shareface = x
									objinterfelems[-1].grain = ynodegrains
									indexobjinterfelems.append([int(x),int(y)])
		self.interfelems = objinterfelems		
		return 1

	# create faces between interface elements
	def CreateFaces(self):
		# they are already arranged in pairs with common face
		for nelem in range(0,len(self.interfelems),2):
			tempelem1 = self.interfelems[nelem]
			tempelem2 = self.interfelems[nelem+1]
			tempnodes1 = tempelem1.nodes
			tempnodes2 = tempelem2.nodes
			# check positions of the common interface nodes in the elements
			posintnodes1 = np.zeros(shape=(4))
			posintnodes2 = np.zeros(shape=(4))
			count = 0
			for n in range(0,8):
				if tempnodes1[n] in tempnodes2:
					posintnodes1[count] = n+1 # index from 1 to 8
					dove = np.where(tempnodes2 == tempnodes1[n])
					tempelem1.corrnodes.append([n+1,int(dove[0])+1]) # index from 1 to 8
					count = count + 1
			count = 0
			for n in range(0,8):
				if tempnodes2[n] in tempnodes1:
					posintnodes2[count] = n+1 # index from 1 to 8
					dove = np.where(tempnodes1 == tempnodes2[n])
					tempelem2.corrnodes.append([n+1,int(dove[0])+1]) # index from 1 to 8
					count = count + 1
			# assign face type to interface element
			if 1 in posintnodes1: # can be type 1, 4 or 5
				if 2 in posintnodes1: # can be type 1 or 5
					if 3 in posintnodes1: # can be type 1 only
						self.interfelems[nelem].facetype = 1
					else: # can be type 5 only
						self.interfelems[nelem].facetype = 5
				else: # can be type 4 only
					self.interfelems[nelem].facetype = 4				
			else: # can be type 2, 3 or 6
				if 2 in posintnodes1: # can be type 3 only
					self.interfelems[nelem].facetype = 3
				else: # can be type 2 or 6
					if 3 in posintnodes1: # can be type 6 only
						self.interfelems[nelem].facetype = 6 
					else: # can be type 2 only
						self.interfelems[nelem].facetype = 2
			# do the same for the pair element
			if 1 in posintnodes2: # can be type 1, 4 or 5
				if 2 in posintnodes2: # can be type 1 or 5
					if 3 in posintnodes2: # can be type 1 only
						self.interfelems[nelem+1].facetype = 1
					else: # can be type 5 only
						self.interfelems[nelem+1].facetype = 5
				else: # can be type 4 only
					self.interfelems[nelem+1].facetype = 4				
			else: # can be type 2, 3 or 6
				if 2 in posintnodes2: # can be type 3 only
					self.interfelems[nelem+1].facetype = 3
				else: # can be type 2 or 6
					if 3 in posintnodes2: # can be type 6 only
						self.interfelems[nelem+1].facetype = 6 
					else: # can be type 2 only
						self.interfelems[nelem+1].facetype = 2
		return 1

	# order cohesive nodes counterclockwise
	def CreateOrderCohNodes(self):
		for n in range(0,len(self.interfelems)):
			self.interfelems[n].OrderCohNodes()
		return 1
									 










