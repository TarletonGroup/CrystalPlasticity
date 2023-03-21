Nicolo Grilli
University of Oxford
AWE project 2020
23 Giugno 2020

Code to generate cohesive elements
at the grain boundaries in 3D

Copy 'Job-1.inp' input file in this folder
Elset must be named 'GRAIN1', 'GRAIN2', et cetera

Run
python3.6 main.py

Old bulk element list from abaqus input file will be in:
'Job-1-bulk-elems.inp'

Old node list from abaqus input file will be in:
'Job-1-node.inp'

Old element sets with the generate keyword will be in:
'Job-1-elset.inp'

New node list will be in file:
'Job-1-node-new.inp'

New bulk element list will be in file:
'Job-1-elems-new.inp'

Interface (zero thickness) element list will be in file:
'Job-1-int-elems.inp'

Maximum triple junctions are handled
A node cannot belong to more than 3 grains
Only C3D8 elements are handled
It is fundamental that nodes at the grain boundary
belongs to bulk elements in two/three grains
and are not already duplicated 

