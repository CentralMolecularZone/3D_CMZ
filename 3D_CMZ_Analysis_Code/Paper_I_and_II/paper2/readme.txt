This directory contains files used to make the figures for 3-D CMZ Paper 2 as well as much of the analysis code. These are not heavily commented or intended to be completely generalized for other data / situations. They are here for transparency of methodology and for adaptation if you find any aspect of it useful. If you find it useful, please cite the paper in any work that utilizes this code or the data products produced in the paper.

The figure files are labelled according to which figure in the text they correspond. In some cases, the figures were then arranged together in a graphical interface to appear as they do in the manuscript. However, these files should contain the complete code for the underlying figures.

The analysis done in paper 2 is also included here in the form of several python scripts. Here was the general workflow for running these scripts:
  1. Run: dendrogram_catalog.py  -- this produces the following megatable elements: 'major_sigma', 'minor_sigma', 'radius', 'area_ellipse', 'area_exact',
          'position_angle', 'x_cen', 'y_cen', 'average_column',
          'median_column', 'peak_column', 'mass' 
  2. Run: add_temperature_to_catalog.py -- this produces the following megatable elements: 'peak_tem', 'mean_tem', 'median_tem'
  3. Run: extract_hcn.py -- this produces the following megatable elements: HCN, HC3N, and HNCO mom0, mom1, and fwhm (as well as producing figures and CO analysis)
  4. Run: add_luminosity_to_catalog.py - this produces the following megatable elements: 'peak_lum', 'mean_lum', 'total_lum', 'total_lum_SFR',
               'peak_lum_no70', 'mean_lum_no70', 'total_lum_no70', 'total_lum_no70_SFR'
  5. Run: MLE_fitting_all_2022update.ipynb -- which does all of the MLE power-law fitting. These are not automatically added to the table
  6. Run: add_free-fall_SFRs_and_compare_with_cmzoom.py -- which adds the free-fall and CMZoom SFRs to the megatable
  7. The mega-table was then configured by hand for latex presentation (colloquiual names etc. were complicated) and the power-law fits were added manually as well.
