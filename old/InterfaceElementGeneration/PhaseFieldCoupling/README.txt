The coupling between phase field in the bulk elements
and fracture properties in the cohesive interface elements
requires to build a look-up table that indicates the bulk elements
that are attached to one cohesive interface elements

This is done by the script:
CoheleBulkMap.py

At the beginning of the script,
set the files with nodes, bulk elements
and interface elements (after reorientation), for instance:

nodefilename = 'Square20umEl0p25-node.inp'
bulkelemfilename = 'Square20umEl0p25-bulk-elems.inp'
intelemfilename = 'Square20umEl0p25-int-elems-new.inp'

The results will be in the file:
CoheleBulkMap.f

This file must be copied to the folder with the umat code
It is already included in the UEL.for with "include" command

WARNING: the more interface elements, the more the compiler
will be slow

