Paper 4 readme

This directory contains files used to make the figures for 3-D CMZ Paper 4 as well as much of the analysis code. These are not heavily commented or intended to be completely generalized for other data / situations. They are here for transparency of methodology and for adaptation if you find any aspect of it useful. If you find it useful, please cite the paper in any work that utilizes this code or the data products produced in the paper.

For any questions or clarifications, please contact Dani Lipman (dani.lipman@uconn.edu)

In the main Paper_IV directory are the main notebooks and python scripts used for the dust extinction analysis.

CMZ_FFT_smoothing.py: performs masking and smoothing of the Spitzer 8um and Herschel 70um maps 

70um_catalog.ipynb: masks out the herschel 70um point source catalog from the 70um map.

flux_fractions.py: performs flux diff and flux ratio method calculations 

8um_extinction_col_den.ipynb: Performs Spitzer 8um extinction column density calculations and creates plots similar to Fig C1

70um_extinction_col_den.ipynb: Performs Herscehl 70um extinction column density calculations and creates plots similar to Fig C1

Sub_masks_ExtN_Calcs.ipynb: Applies the cloud submasks from Paper III and recalculates all extinction methods for the sub masks.


The subfolders contain code for further analysis and main figures.
3D_CMZ_Analysis_Code/Paper_IV/paper_figs/ contains code that was used to create the main overview figures in paper 4, including Figs 1, 3 and 4. 

3D_CMZ_Analysis_Code/Paper_IV/LBV_plots/ contains code and files used for the lbv plots and comparisons between methods (Figs 7, 8, 9, 10, 11).

3D_CMZ_Analysis_Code/Paper_IV/4d_comparison contains code used for the percent agreement method (Table 2 and Fig 12)

