# Codes_DOL_Mechanisms_Nonclonal
This repository is for open access to source codes and data for the manuscript **The evolution of mechanisms to divide labour in non-clonal population**

In **analytical_scripts** are the files for generating Fig.3 and Fig.S1-3. The enclosed scripts do not label axes or results. For this information please refer to the figures as they occur in the main text and supplementary information. Please include WK.m in the same folder in order to run: FiguresS1BandC.m, FigureS2C.m and FigureS2D.m.

In **simulations** are the files for Fig.S4: */codes* contains source codes for running time series and the heatmap, while */data* contains raw data and plotting codes. In terms of the source codes, please *do not* change the relative location of these files as the random number generators files (dSFMT) are required for running, and **please read** the notes at the beginning of each file, containing the key parameters and notes. For the data files, the two heatmaps of coordination levels and proportion of helpers are plotted from *level_sumary.csv*, you may use the files in *R_files* to replot the figures. If you wish to generate your own *level_summary.csv*, please assign the folders as in *data/Fig4_levelplot* and run *cord_ss.c* for each group size separately. The *merge_summaries.R* would help you create the *level_summary.csv* file. 

Please email *ming.liu@zoo.ox.ac.uk* if there is any problem, thanks! (Ming Liu)
