# Codes_DOL_Mechanisms_Nonclonal
This repository is for open access to source codes and data for the manuscript **The evolution of mechanisms to divide labour in non-clonal population**

**Codes** contains source codes of all models: *analytical_q_residents* is the main analytical model, *analytical_q_optimal* is the analytical model in section 2 of supplementary information, and *simulations* are the codes for running the simulative model. For */simulations* source codes, please *make sure* the random number generators files (dSFMT) are correctly located in the *.c* files because those files are required for executing the program. Also, **please read** the notes at the beginning of each file as they provide descriptions of key parameters and instructions for running. 

In **sims_data** are the data used for plotting simulation results, including supplementary figures. If you wish to generate your own *level_summary.csv* for the level plot, please run *hsta_ns.c* in *codes/simulations/repeats* for each group size separately with the same number of repetition. *Merge_summaries.R* would help you merge data to export *level_summary.csv*. 

Please email *ming.liu@zoo.ox.ac.uk* if there is any problem, thanks! (Ming Liu)
