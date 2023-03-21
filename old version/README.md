UMAT for Abaqus written by Nicolò Grilli and Ed Tarleton based on UEL by Fionn Dunne.

Grilli, N., Tarleton, E. & Cocks, A.C.F. Coupling a discrete twin model with cohesive elements to understand twin-induced fracture. Int J Fract 227, 173–192 (2021). https://doi.org/10.1007/s10704-020-00504-9

# Usage instructions:

define a user material in Abaqus with
125 state dependent variables (SDV)
11 material constants

The first constant indicates the crystal type:

0 = HCP
1 = BCC
2 = BCC
3 = Carbide
4 = Olivine
5 = Orthorombic

The constants 2-10 contains the components of the rotation matrix
that transforms a vector from the crystal reference frame
to the sample (Abaqus) reference frame
The order of the components in the input file must be
R11, R12, R13, R21, R22, R23, R31, R32, R33

The 11th constant is the grain index
different for different grains
Therefore a polycrystal should be made with
different materials in abaqus

The parameters that must be set are in the variable declaration part
of the following files:
umat.for
kmat.f
kMaterialParam.f

The total number of bulk elements and cohesive interface elements 
must be set in:
mycommon.f

If the twins are activates, it is necessary to set the number
of neighbouring points (NUpDown) in mycommon.f
and then set ArrayNUpDown as the product:
NUpDown * nElements * nintpts
A preliminary simulation must be run,
if the number of neighbouring points is not high enough
for the specific mesh, warning will be given in the .log file
The last warning line will contain the number that must be
assigned to NUpDown

IMPORTANT: when using discrete twin model
be sure no more such warnings are present in the .log
for reliable simulation results

kMaterialParam.f includes the elastic and plastic parameters
for different materials. If a new material is needed,
new parameters and material name must be introduced in this file

kRhoTwinInit.f is preset, but can be modified
if it is necessary to introduce a specific initial value
(for instance space dependent) of the dislocation density 
or twin phase field





