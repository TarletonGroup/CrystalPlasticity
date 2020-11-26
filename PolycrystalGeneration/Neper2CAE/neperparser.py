# Nicolo Grilli
# University of Oxford
# AWE project 2020
# 3 luglio 2020

# class representing neper tess file
# and functions for parsing

class NeperParser:

	# prefissofile is the prefix of neper .tess file
	def __init__(self,prefissofile):
		self.prefissofile = prefissofile
		self.filename = prefissofile + '.tess'

	# read Euler angles from file
	def ReadEulerAngles(self):
		fid = open(self.filename,'r')
		euleranglesfilename = self.prefissofile + '-eulerangles.tess'
		fout = open(euleranglesfilename,'w')

		flagcheckekeyword = False # e keyword after *ori
		flagprintnextline = False # next line should be printed or not

		for line in fid:
			if (flagprintnextline):
				if '*' in line: # end of Euler angles
					flagprintnextline = False
				else:
					fout.write(line)
			if (flagcheckekeyword):
				flagprintnextline = True
				flagcheckekeyword = False
			if '*ori' in line: # this is the Euler angles keyword
				flagcheckekeyword = True

		fid.close()
		fout.close()

		return 1

	# read vertices that will become datum points for abaqus
	def ReadVertices(self):
		fid = open(self.filename,'r')
		datumpointsfilename = self.prefissofile + '-datumpoints.tess'
		fout = open(datumpointsfilename,'w')

		flagprintnextline = False # next line should be printed or not
		flagNvertices = False # number of vertices is reported after vertices keyword

		for line in fid:
			if (flagprintnextline):
				if '*' in line: # end of vertices
					flagprintnextline = False
				else:
					fout.write(line)
			if (flagNvertices):
				nverttemp = line.split()
				print('Number of vertices in tess file:' + '\n')
				print(int(nverttemp[0]))
				flagNvertices = False
				flagprintnextline = True
			if '**vertex' in line: # this is the vertices keyword
				flagNvertices = True

		fid.close()
		fout.close()

		return 1

	# read edges that will be used by abaqus to partition face 
	def ReadEdges(self):
		fid = open(self.filename,'r')
		edgesfilename = self.prefissofile + '-edges.tess'
		fout = open(edgesfilename,'w')

		flagprintnextline = False # next line should be printed or not
		flagNedges = False # number of edges is reported after edge keyword

		for line in fid:
			if (flagprintnextline):
				if '*' in line: # end of edges
					flagprintnextline = False
				else:
					fout.write(line)
			if (flagNedges):
				nedgetemp = line.split()
				print('Number of edges in tess file:' + '\n')
				print(int(nedgetemp[0]))
				flagNedges = False
				flagprintnextline = True
			if '**edge' in line:
				flagNedges = True

		fid.close()
		fout.close()

		return 1















