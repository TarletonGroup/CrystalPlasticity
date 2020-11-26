# Nicolo Grilli
# University of Oxford
# AWE project 2020
# 5 giugno 2020

# class representing a COH3D8 element
# at the boundary between grains

import numpy as np
from interfelement import InterfElement

class CohElement:

	def __init__(self,elem1,elem2):
		# a cohesive elements is at the interface
		# of two interface elements
		self.elem1 = elem1
		self.elem2 = elem2
		# nodes
		self.nodes = np.zeros(shape=(8))
		# index must follow the indices of bulk elements
		self.index = 0

	# assign nodes based on the input interface elements
	def AssignNodes(self):
		# follow the order of the first interface element
		# so first interface element face will be counterclockwise
		# and second interface element face will be clockwise
		temp1order = np.int_(self.elem1.cohnodeorder)
		temp1corr = self.elem1.corrnodes
		tempnodes = np.zeros(shape=(8))
		for n in range(0,4):
			tempnodes[n] = int(self.elem1.nodes[temp1order[n]-1])
			for search in range(0,4):
				if (temp1corr[search][0] == temp1order[n]):
					tempindex2 = int(temp1corr[search][1])
					tempnodes[n+4] = int(self.elem2.nodes[tempindex2-1])
		self.nodes = np.int_(tempnodes)
		return 1


	






