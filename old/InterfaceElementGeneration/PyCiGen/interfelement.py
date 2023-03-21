# Nicolo Grilli
# University of Oxford
# AWE project 2020
# 5 giugno 2020

# class representing a C3D8 element
# with four nodes on the grain boundary

import numpy as np
from element import Element

class InterfElement(Element):

	def __init__(self,index,nodes):
		super().__init__(index,nodes)
		# which grain does this element belong?
		self.grain = 0
		# index of the other element that shares
		# a face on the grain boundary with self
		self.shareface = 0
		# 6 possible grain boundary faces on C3D8 elements
		# 1 -> 1,2,3,4
		# 2 -> 5,6,7,8
		# 3 -> 2,3,6,7
		# 4 -> 1,4,5,8
		# 5 -> 1,2,5,6
		# 6 -> 3,4,7,8
		self.facetype = 0
		# store correspondence of nodes
		# this is a list of [n1,n2] in which n1 is the
		# index of the node in this elements and n2 is the index
		# of the corresponding node in the other elements sharing the face
		# indices go from 1 to 8
		self.corrnodes = []
		# order of the nodes in the cohesive
		# element on this face
		self.cohnodeorder = np.zeros(shape=(4))

	# define order of the nodes on the cohesive element
	# self.facetype must be assigned
	def OrderCohNodes(self):
		if (self.facetype > 0):
			if (self.facetype == 1):
				self.cohnodeorder = np.array([4,3,2,1])
			if (self.facetype == 2):
				self.cohnodeorder = np.array([5,6,7,8])
			if (self.facetype == 3):
				self.cohnodeorder = np.array([2,3,7,6])
			if (self.facetype == 4):
				self.cohnodeorder = np.array([1,5,8,4])
			if (self.facetype == 5):
				self.cohnodeorder = np.array([1,2,6,5])
			if (self.facetype == 6):
				self.cohnodeorder = np.array([3,4,8,7])
		else:
			print('error: facetype not assigned')
		return 1


	






