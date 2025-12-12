This is the README file for Paper V. Below are descriptions of the code necessary for the analysis and figures presented in Paper V.

Note, the code has been lightly commented for transparency. For any questions or issues, please contact Dani R. Lipman.



--- 

coords.py is a utility file with functions needed to project coordinates from cartesian to galactic frames.

star_counts_modeling.ipynb creates star count background model map shown in the third panel of Fig 2. 

star_count_NF_method.ipynb uses the model map for the stat count ratio method, and produces Figure 2, and starcounts_tab.tex.

make_PPDFs.ipynb produces the prior and posterior PPDFs seen in Fig Fig 3 and A.1. It utilizes values from synth_table.tex, starcounts_tab.tex, nogueras_measures.tex, and xray_methods.tex to create priors. The values of the posterior distributions are reported in Table 2.

x2_fitting_mahalanobis_wasserstein.ipynb performs the x2 orbital parameter grid fitting and ptop-down projection/distance estimates. Figures 5 and 6 are created here.

fitting_comparisons.ipynb creates Figures 4, 7, and 8. 

---

