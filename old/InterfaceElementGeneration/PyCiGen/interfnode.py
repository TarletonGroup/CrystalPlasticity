# Nicolo Grilli
# University of Oxford
# AWE project 2020
# 4 giugno 2020

# class representing a node 
# on a grain boundary

import numpy as np
from node import Node

class Interfnode(Node):

	def __init__(self,index,coords):
		super().__init__(index, coords)
		# to how many grains does this node belong?
		# can be 2 or 3
		self.multiplicity = 0
		# indices of the twin/triplet duplicated nodes
		self.twinindex = 0
		self.tripletindex = 0
		# after duplication this node will belong to a single grain
		self.grain = 0
		# after duplication connectivity will change
		# only elements in the same grain will be present
		self.newconnectivity = np.zeros(shape=(0))

