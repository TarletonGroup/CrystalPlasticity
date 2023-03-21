# Nicolo Grilli
# University of Oxford
# AWE project 2020
# 4 giugno 2020

# class representing a twin (or triplet) node
# index will be = number of node + index of twins
# it will be the index in the new input file
# triplets will be duplicated twice
# and will be next to each other in the new input file 

import numpy as np
from node import Node

class TwinNode(Node):

	def __init__(self,index,coords):
		super().__init__(index, coords)
		# to how many grains does this node belong?
		# can be 2 or 3
		self.multiplicity = 0
		self.originalnode = 0 # original abaqus index of the node
		# this node will belong to a single grain
		self.grain = 0
		# after duplication connectivity will change
		# only elements in the same grain will be present
		self.newconnectivity = np.zeros(shape=(0))

