import runpy

from paths import rootpath

runpy.run_path('dendrogram_catalog.py')
runpy.run_path('add_temperature_to_catalog.py')
runpy.run_path('extract_hcn.py')
runpy.run_path('add_luminosity_to_catalog.py')
#runpy.run_path('add_powerlaws_to_table.py')
#runpy.run_path('make_tex_file.py')
#runpy.run_path('make_short_table.py')
