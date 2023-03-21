Nicolo Grilli
University of Oxford
AWE project 2020
1 Agosto 2020

Code to link Neper polycrystal generator
with Abaqus CAE

In main.py, set the parameters:
width
height
depth
in arbitrary units

These will be the dimensions of the parallelepiped
geometry along the x,y,z axes generated in Abaqus

In main.py, set the parameter:
Ngrains

This will be the number of grains

Execute:

python3.6 main.py

Abaqus CAE will be opened and columnar grains generated
and materials and sections assigned to them

Grain growth model is used by default
but custom Neper options can be added in main.py

Custom material constants for the grains 
can be introduced in GeneratePolycrystal.py

Requires Neper and Abaqus installed,
otherwise the generated files can be moved to another
machine with Abaqus installed and the command

abaqus cae script=GeneratePolycrystal.py

can be exectuted independently on that machine in the same folder

