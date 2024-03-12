Eralp Demir
March 10th, 2024

This example contains the input parameters sections 3.6 on irradiation effects.

The input file contains only the analysis for: 
- irradiation model-1 (simplified model)
- irradiation model-2 (rigorous model) - material type HCP


The userinputs.f file needs to be copied to the UMAT directory that contains Fortan files to ensure the selected input parameters are used.
For example (some important ones):
- numel = 8
- eltyp = 'C3D8'
- maxnslip = 30

