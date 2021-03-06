# Codes_DOL_Mechanisms_Nonclonal
This repository is for open access to source codes and data for the manuscript entitled **"The evolution of mechanisms to divide labour in non-clonal populations"**.

**Codes** contains source codes of all models: *analytical/q_residents* is the main analytical model, *analytical/q_optimal* is the analytical model in section 2 of supplementary information, *analytical/var_theta* is the codes for Fig S1, and *simulations/* are the codes for running the simulative model. 

For the C files in *codes/simulations*, **please** make sure the location of random number generators files match the description in C files (i.e., *#include "../../dSFMT-src-2.2.3/dSFMT.c"*) because those files are essential for running. Also, **please** read the notes at the beginning of each file as they provide descriptions of key parameters and instructions for running. 

In **sims_data** are the data used for plotting simulation results, including supplementary figures. If you wish to generate your own *level_summary.csv* for the level plot, please run *hsta_ns.c* in *codes/simulations/repeats* for each group size separately with the same number of repetition. *Merge_summaries.R* would help you merge data to export *level_summary.csv*. 

Please email *ming.liu@zoo.ox.ac.uk* if there is any problem, thanks! (Ming Liu)
