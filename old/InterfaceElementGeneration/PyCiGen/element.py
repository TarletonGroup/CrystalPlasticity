# Nicolo Grilli
# University of Oxford
# AWE project 2020
# 4 giugno 2020

# class representing a C3D8 element

import numpy as np

class Element:

	def __init__(self,index,nodes):
		self.index = int(index)
		self.nodes = np.int_(nodes)
		# which grain does this element belong?
		self.grain = 0


	

