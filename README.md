# Analysis of STEP UP with functional data analysis

Code accompanying "Quantifying Physical Activity Intervention Effects via Functional Regression".

Public data is not available.

This repository contains the following code files:
- fpca.R: Conduct function principal components analysis (FPCA) and run a regression of eigenfunction scores on covariates
- fosr_fpcareg.R: Function-on-a-scalar regression, FPCA + regression as two-step model, inference for both methods via bootstrap
- fofr.R: Function-on-function regression and inference via bootstrap
- simulation_study.R: Simulation study of FoSR in the intervention period