# Nicolo Grilli
# University of Oxford
# AWE project 2020
# 4 giugno 2020

# class representing a node of a C3D8 element

import numpy as np

class Node:

	def __init__(self,index,coords):
		self.index = int(index)
		self.coords = coords
		# indices of the elements that contain this node (abaqus index)
		self.connectivity = np.zeros(shape=(0))
