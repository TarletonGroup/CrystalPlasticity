Eralp Demir
March 10th, 2024

This example contains the input parameters sections 3.6 on irradiation effects.

The input files in two folders contain only the analysis for: 
- irradiation model-1 (simplified model)
- irradiation model-2 (rigorous model) - material type HCP


The userinputs.f file needs to be modified:

- To model HCP materials set "maxnslip = 30"

