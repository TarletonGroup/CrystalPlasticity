# Nicolo Grilli
# University of Oxford
# AWE project 2020
# 5 giugno 2020

# class representing the mesh
# with interface elements

import numpy as np
from numpy import genfromtxt
from mesh import Mesh
from element import Element
from node import Node
from interfnode import Interfnode
from twinnode import TwinNode
from cohelement import CohElement

class NewMesh:

	# construct based on the old mesh
	def __init__(self,mesh):
		self.nodes = mesh.nodes
		self.bulkelements = mesh.bulkelements
		self.interfnodes = mesh.interfnodes
		self.interfelems = mesh.interfelems
		self.twinnodes = []
		# initialize new elements as the old ones, then modify
		self.newbulkelems = mesh.bulkelements
		# initialize new interface elements as the old ones, then modify
		self.newinterfelems = mesh.interfelems
		# cohesive elements
		self.cohelems = []

	# duplicate nodes and create twin/triplet nodes
	def CreateTwinNodes(self):
		objtwinnodes = []
		Nnodes = len(self.nodes) # original number of nodes
		count = Nnodes # increasing index for the duplicated nodes
		for objinterfnode in self.interfnodes:
			if (objinterfnode.multiplicity == 2): # this is a twin
				count = count + 1
				objtwinnodes.append(TwinNode(count,objinterfnode.coords))
				objtwinnodes[-1].originalnode = int(objinterfnode.index)
				objinterfnode.twinindex = count
				objtwinnodes[-1].connectivity = objinterfnode.connectivity
				objtwinnodes[-1].multiplicity = objinterfnode.multiplicity
				for elemindex in objtwinnodes[-1].connectivity:
					if (self.bulkelements[elemindex-1].grain != objinterfnode.grain):
						objtwinnodes[-1].grain = self.bulkelements[elemindex-1].grain
						break
				# assign new connectivity of twin node
				for elemindex in objtwinnodes[-1].connectivity:
					if (self.bulkelements[elemindex-1].grain == objtwinnodes[-1].grain):
						objtwinnodes[-1].newconnectivity = np.int_(np.append(objtwinnodes[-1].newconnectivity,elemindex))
			if (objinterfnode.multiplicity == 3): # this is a triplet
				count = count + 1
				objtwinnodes.append(TwinNode(count,objinterfnode.coords))
				objtwinnodes[-1].originalnode = int(objinterfnode.index)
				objinterfnode.twinindex = count
				objtwinnodes[-1].connectivity = objinterfnode.connectivity
				objtwinnodes[-1].multiplicity = objinterfnode.multiplicity
				for elemindex in objtwinnodes[-1].connectivity:
					if (self.bulkelements[elemindex-1].grain != objinterfnode.grain):
						objtwinnodes[-1].grain = self.bulkelements[elemindex-1].grain
						break
				# assign new connectivity of twin node
				for elemindex in objtwinnodes[-1].connectivity:
					if (self.bulkelements[elemindex-1].grain == objtwinnodes[-1].grain):
						objtwinnodes[-1].newconnectivity = np.int_(np.append(objtwinnodes[-1].newconnectivity,elemindex))
				count = count + 1
				objtwinnodes.append(TwinNode(count,objinterfnode.coords))
				objtwinnodes[-1].originalnode = int(objinterfnode.index)
				objinterfnode.tripletindex = count
				objtwinnodes[-1].connectivity = objinterfnode.connectivity
				objtwinnodes[-1].multiplicity = objinterfnode.multiplicity
				for elemindex in objtwinnodes[-1].connectivity:
					if (self.bulkelements[elemindex-1].grain != objinterfnode.grain):
						if (self.bulkelements[elemindex-1].grain != objtwinnodes[-2].grain):
							objtwinnodes[-1].grain = self.bulkelements[elemindex-1].grain
							break
				# assign new connectivity of triplet node
				for elemindex in objtwinnodes[-1].connectivity:
					if (self.bulkelements[elemindex-1].grain == objtwinnodes[-1].grain):
						objtwinnodes[-1].newconnectivity = np.int_(np.append(objtwinnodes[-1].newconnectivity,elemindex))
		self.twinnodes = objtwinnodes
		return 1

	# assign duplicated nodes to elements
	def MakeNewBulkElems(self):
		for node in self.twinnodes:
			tempconnectivity = node.newconnectivity
			# only elements connected to twin nodes must be changed
			for elem in tempconnectivity: # elem and node belong to same grain already
				count = 0
				for tempnode in self.newbulkelems[elem-1].nodes:
					# if this node is an original node, then substitute with twin node
					if (self.newbulkelems[elem-1].nodes[count] == node.originalnode):
						self.newbulkelems[elem-1].nodes[count] = node.index
					count = count + 1	
		return 1

	# assign duplicated nodes to interface elements
	# just copy from new bulk elements modified earlier
	def MakeNewInterfElems(self):
		for nelem in range(0,len(self.newinterfelems)):
			bulkindex = self.newinterfelems[nelem].index # index of the corresponding bulk element
			self.newinterfelems[nelem].nodes = np.int_(self.newbulkelems[bulkindex-1].nodes)
		return 1

	# build cohesive elements
	# self.newinterfelems have the new nodes after the application of MakeNewInterfElems
	# the nodes of new bulk elements must be used
	def CreateCohesiveElems(self):
		objcohel = []
		# they are already arranged in pairs with common face
		# index of the interface elements is the
		# same as the index of the new bulk elements
		count = len(self.newbulkelems)
		for nelem in range(0,len(self.interfelems),2):
			tempelem1 = self.newinterfelems[nelem]
			tempelem2 = self.newinterfelems[nelem+1]
			objcohel.append(CohElement(tempelem1,tempelem2))
			objcohel[-1].AssignNodes()
			count = count + 1
			objcohel[-1].index = count
		self.cohelems = objcohel
		return 1

	# write new elements file
	def WriteNewElemsFile(self,prefissofile):
		objelems = self.newbulkelems
		Nelems = len(self.newbulkelems) # number of elements
		newelems = np.zeros(shape=(Nelems,9))
		for n in range(0,Nelems):
			newelems[n,0] = int(objelems[n].index)
			newelems[n,1:9] = np.int_(objelems[n].nodes)
		np.savetxt(prefissofile + '-bulk-elems-new.inp', newelems, fmt='%d,%d,%d,%d,%d,%d,%d,%d,%d')
		return 1

	# write new node file
	def WriteNewNodeFile(self,prefissofile):
		objnodes = self.nodes
		Nnodes = len(self.nodes) # original number of nodes
		objtwinnodes = self.twinnodes
		Ntwinnodes = len(self.twinnodes)
		newnode = np.zeros(shape=(Nnodes+Ntwinnodes,4))
		for n in range(0,Nnodes):
			tempnode = objnodes[n]
			newnode[n,0] = tempnode.index
			newnode[n,1:4] = tempnode.coords
		for n in range(0,Ntwinnodes):
			tempnode = objtwinnodes[n]
			newnode[n+Nnodes,0] = tempnode.index
			newnode[n+Nnodes,1:4] = tempnode.coords
		np.savetxt(prefissofile + '-node-new.inp', newnode, fmt='%d,%f,%f,%f')

	# write cohesive elements
	def WriteCohElems(self,prefissofile):
		objcohelems = self.cohelems
		Ncohelems = len(self.cohelems) # number of cohesive elements
		cohelemarray = np.zeros(shape=(Ncohelems,9))
		for n in range(0,Ncohelems):
			tempcohelem = objcohelems[n]
			cohelemarray[n,0] = int(tempcohelem.index)
			cohelemarray[n,1:9] = np.int_(tempcohelem.nodes)
		np.savetxt(prefissofile + '-int-elems.inp', cohelemarray, fmt='%d,%d,%d,%d,%d,%d,%d,%d,%d')
		return 1

	# write dummy elements to visualise UEL state variables
	def WriteDummyElems(self,prefissofile):
		# they are already arranged in pairs with common face
		# index of the dummy element runs first over the first neighbouring element
		# of each interface elements, then over the second neighbouring element
		# of each interface element
		Ncohelems = len(self.cohelems)
		count = len(self.newbulkelems) + Ncohelems
		count2 = 0
		# same number of dummy elements as interface elements
		Ndummyelems = 2*Ncohelems
		dummyelemsarray = np.zeros(shape=(Ndummyelems,9))
		for nelem in range(0,Ndummyelems,2):
			tempelem1 = self.newinterfelems[nelem]
			tempelem2 = self.newinterfelems[nelem+1]
			dummyelemsarray[count2,1:9] = np.int_(tempelem1.nodes)
			dummyelemsarray[count2+Ncohelems,1:9] = np.int_(tempelem2.nodes)
			count = count + 1
			dummyelemsarray[count2,0] = count
			dummyelemsarray[count2+Ncohelems,0] = count+Ncohelems
			count2 = count2 + 1
		np.savetxt(prefissofile + '-dummy-elems.inp', dummyelemsarray, fmt='%d,%d,%d,%d,%d,%d,%d,%d,%d')
		return 1
		 











