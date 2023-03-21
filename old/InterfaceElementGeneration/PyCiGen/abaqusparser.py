# Nicolo Grilli
# University of Oxford
# AWE project 2020
# 23 Giugno 2020

# class representing abaqus input file
# and functions for parsing

class AbaqusParser:

	# prefissofile is the prefix of abaqus .inp file
	def __init__(self,prefissofile):
		self.prefissofile = prefissofile
		self.filename = prefissofile + '.inp'

	# read element sets from file
	# only generate keyword is accepted
	def ReadElsets(self):
		fid = open(self.filename,'r')
		elsetsfilename = self.prefissofile + '-elsets.inp'
		fout = open(elsetsfilename,'w')

		flagprintnextline = False # next line should be printed or not

		for line in fid:
			if (flagprintnextline):
				fout.write(line)
				flagprintnextline = False
			if (line[0:6] == '*Elset'): # this is an elset
				posgrain = line.find('GRAIN') # -1 if this is not a grain
				posgenerate = line.find('generate') # generate is only format accepted
				if (posgrain >= 0): # this is a grain
					if (posgenerate >= 0):
						fout.write(line)
						flagprintnextline = True # print next line containing the elements start/end
					else:
						print('generate keyword is not present' + '\n')
						print('generate keyword is the only accepted format' + '\n')
						exit()

		fid.close()
		fout.close()

		return 1

	# read nodes from file
	def ReadNodes(self):
		fid = open(self.filename,'r')
		nodesfilename = self.prefissofile + '-node.inp'
		fout = open(nodesfilename,'w')

		flagfoundnodes = False # block with nodes has been found

		for line in fid:
			if (flagfoundnodes): # this line may contain nodes
				if (line[0] == '*'): # reached the end of nodes
					flagfoundnodes = False
					break # needed because you may find also *Node Output
				else:
					fout.write(line)
			if (line[0:5] == '*Node'): # these are the nodes
				flagfoundnodes = True

		fid.close()
		fout.close()

		return 1

	# read bulk elements from file
	def ReadBulkElems(self):
		fid = open(self.filename,'r')
		elemsfilename = self.prefissofile + '-bulk-elems.inp'
		fout = open(elemsfilename,'w')

		flagfoundelems = False # block with elements has been found

		for line in fid:
			if (flagfoundelems): # this line may contain elements
				if (line[0] == '*'): # reached the end of elements
					flagfoundelems = False
					break # needed because you may find also *Element Output
				else:
					fout.write(line)
			if (line[0:8] == '*Element'): # these are the elements
				flagfoundelems = True

		fid.close()
		fout.close()

		return 1







