# Nicolo Grilli
# University of Oxford
# AWE project 2020
# 26 Giugno 2020

# generate sets where boundary conditions are applied
# based on the node coordinates

import numpy as np
from numpy import genfromtxt

# read node file
nodefilename = 'Square20umEl0p25-node.inp'
nodefile = genfromtxt(nodefilename,dtype=float,delimiter=",",skip_header=0,skip_footer=0)

Xmin = None
Xmax = None
Ymin = None
Ymax = None
Zmin = None

# search for min/max values of the coordinates
for node in range(0,len(nodefile)):
	if Xmin is None:
		Xmin = nodefile[node,1]
	else:
		Xmin = min(Xmin,nodefile[node,1])
	if Xmax is None:
		Xmax = nodefile[node,1]
	else:
		Xmax = max(Xmax,nodefile[node,1])
	if Zmin is None:
		Zmin = nodefile[node,3]
	else:
		Zmin = min(Zmin,nodefile[node,3])
	if Ymin is None:
		Ymin = nodefile[node,2]
	else:
		Ymin = min(Ymin,nodefile[node,2])
	if Ymax is None:
		Ymax = nodefile[node,2]
	else:
		Ymax = max(Ymax,nodefile[node,2])

print('Minimum value of X coordinate is:')
print(Xmin)
print('Maximum value of X coordinate is:')
print(Xmax)
print('Minimum value of Y coordinate is:')
print(Ymin)
print('Maximum value of Y coordinate is:')
print(Ymax)
print('Minimum value of Z coordinate is:')
print(Zmin)

# build node sets for the boundary
nodesOrigin = np.array([],dtype='uint64')
nodesXmin = np.array([],dtype='uint64')
nodesXmax = np.array([],dtype='uint64')
nodesYmin = np.array([],dtype='uint64')
nodesYmax = np.array([],dtype='uint64')
nodesZmin = np.array([],dtype='uint64')

# build boundary sets
for node in range(0,len(nodefile)):
	if (nodefile[node,1] == Xmin and nodefile[node,2] == Ymin and nodefile[node,3] == Zmin): # Origin
		nodesOrigin = np.append(nodesOrigin,np.uint64(nodefile[node,0]))
	if (nodefile[node,1] == Xmin): # Xmin
		nodesXmin = np.append(nodesXmin,np.uint64(nodefile[node,0]))
	if (nodefile[node,1] == Xmax): # Xmax
		nodesXmax = np.append(nodesXmax,np.uint64(nodefile[node,0]))
	if (nodefile[node,2] == Ymin): # Ymin
		nodesYmin = np.append(nodesYmin,np.uint64(nodefile[node,0]))
	if (nodefile[node,2] == Ymax): # Ymax
		nodesYmax = np.append(nodesYmax,np.uint64(nodefile[node,0]))
	if (nodefile[node,3] == Zmin): # Zmin
		nodesZmin = np.append(nodesZmin,np.uint64(nodefile[node,0]))

# write node sets for the assembly
fid = open('Job-1-node-sets.inp','w')

fid.write("** ASSEMBLY" + "\n")
fid.write("**" + "\n")
fid.write("*Assembly, name=Assembly" + "\n")
fid.write("**" + "\n")
fid.write("*Instance, name=Part-1-1, part=Part-1" + "\n")
fid.write("*End Instance" + "\n")
fid.write("**" + "\n")

fid.write("*Nset, nset=ORIGIN, instance=Part-1-1" + "\n")

count = 0

for x in range(0,len(nodesOrigin)): # Origin
	if (count == 8):
		fid.write("\n")
		count = 0
	fid.write("%d," % nodesOrigin[x])
	count = count + 1 

fid.write("\n")
fid.write("*Nset, nset=X0, instance=Part-1-1")
fid.write("\n")

count = 0

for x in range(0,len(nodesXmin)): # X0
	if (count == 8):
		fid.write("\n")
		count = 0
	fid.write("%d," % nodesXmin[x])
	count = count + 1 

fid.write("\n")
fid.write("*Nset, nset=Xmax, instance=Part-1-1")
fid.write("\n")

count = 0

for x in range(0,len(nodesXmax)): # Xmax
	if (count == 8):
		fid.write("\n")
		count = 0
	fid.write("%d," % nodesXmax[x])
	count = count + 1 

fid.write("\n")
fid.write("*Nset, nset=Y0, instance=Part-1-1")
fid.write("\n")

count = 0

for x in range(0,len(nodesYmin)): # Ymin
	if (count == 8):
		fid.write("\n")
		count = 0
	fid.write("%d," % nodesYmin[x])
	count = count + 1 

fid.write("\n")
fid.write("*Nset, nset=Ymax, instance=Part-1-1")
fid.write("\n")

count = 0

for x in range(0,len(nodesYmax)): # Ymax
	if (count == 8):
		fid.write("\n")
		count = 0
	fid.write("%d," % nodesYmax[x])
	count = count + 1 

fid.write("\n")
fid.write("*Nset, nset=Z0, instance=Part-1-1")
fid.write("\n")

count = 0

for x in range(0,len(nodesZmin)): # Zmin
	if (count == 8):
		fid.write("\n")
		count = 0
	fid.write("%d," % nodesZmin[x])
	count = count + 1 

fid.close()










