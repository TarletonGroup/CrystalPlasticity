Eralp Demir
March 3rd, 2024

This example contains the input parameters sections 3.5 on backstress.

The input file contains only the analysis for the backstres model-1 (A-F model)


The userinputs.f file needs to be copied to the UMAT directory that contains Fortan files to ensure the selected input parameters are used.
For example (some important ones):
- numel = 64000
- eltyp = 'C3D8'
- backstressmodel = 1

