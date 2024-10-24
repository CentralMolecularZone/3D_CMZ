 #!/bin/bash

time nice idl<<EOF

.r convmakesmallimages.pro
.r makebg_gauss_iter.pro
.r gaussfit_iter_sedfit_wrapper_no70beta175_errors.pro


convmakesmallimages, 35.9, 7.8, 'l015', /regridimages, /convimages
convmakesmallimages, 25.2, 5.5, 'l015', /regridimages, /convimages


convmakesmallimages, 35.9, 7.8, 'l017', /regridimages, /convimages
convmakesmallimages, 25.2, 5.5, 'l017', /regridimages, /convimages

convmakesmallimages, 35.9, 7.8, 'l019', /regridimages, /convimages
convmakesmallimages, 25.2, 5.5, 'l019', /regridimages, /convimages


convmakesmallimages, 35.9, 7.8, 'l338', /regridimages, /convimages
convmakesmallimages, 25.2, 5.5, 'l338', /regridimages, /convimages

convmakesmallimages, 35.9, 7.8, 'l341', /regridimages, /convimages
convmakesmallimages, 25.2, 5.5, 'l341', /regridimages, /convimages

convmakesmallimages, 35.9, 7.8, 'l343', /regridimages, /convimages
convmakesmallimages, 25.2, 5.5, 'l343', /regridimages, /convimages

convmakesmallimages, 35.9, 7.8, 'l345', /regridimages, /convimages
convmakesmallimages, 25.2, 5.5, 'l345', /regridimages, /convimages

convmakesmallimages, 35.9, 7.8, 'l347', /regridimages, /convimages
convmakesmallimages, 25.2, 5.5, 'l347', /regridimages, /convimages

convmakesmallimages, 35.9, 7.8, 'l354', /regridimages, /convimages
convmakesmallimages, 25.2, 5.5, 'l354', /regridimages, /convimages

convmakesmallimages, 35.9, 7.8, 'l356', /regridimages, /convimages
convmakesmallimages, 25.2, 5.5, 'l356', /regridimages, /convimages

convmakesmallimages, 35.9, 7.8, 'l358', /regridimages, /convimages
convmakesmallimages, 25.2, 5.5, 'l358', /regridimages, /convimages



makebg_gauss_iter, 4.25, 'l000', /doiter_var, /makelabel, /make36sizebg, /make25sizebg

gaussfit_iter_sedfit_wrapper_no70beta175_errors, 'l000', 'var', /dobackground, /res36, /res25



makebg_gauss_iter, 4.25, 'l002', /doiter_var, /makelabel, /make36sizebg, /make25sizebg

gaussfit_iter_sedfit_wrapper_no70beta175_errors, 'l002', 'var', /dobackground, /res36, /res25



makebg_gauss_iter, 4.25, 'l004', /doiter_var, /makelabel, /make36sizebg, /make25sizebg

gaussfit_iter_sedfit_wrapper_no70beta175_errors, 'l004', 'var', /dobackground, /res36, /res25



makebg_gauss_iter, 4.25, 'l006', /doiter_var, /makelabel, /make36sizebg, /make25sizebg

gaussfit_iter_sedfit_wrapper_no70beta175_errors, 'l006', 'var', /dobackground, /res36, /res25




makebg_gauss_iter, 4.25, 'l008', /doiter_var, /makelabel, /make36sizebg, /make25sizebg

gaussfit_iter_sedfit_wrapper_no70beta175_errors, 'l008', 'var', /dobackground, /res36, /res25




makebg_gauss_iter, 4.25, 'l010', /doiter_var, /makelabel, /make36sizebg, /make25sizebg

gaussfit_iter_sedfit_wrapper_no70beta175_errors, 'l010', 'var', /dobackground, /res36, /res25




makebg_gauss_iter, 4.25, 'l013', /doiter_var, /makelabel, /make36sizebg, /make25sizebg

gaussfit_iter_sedfit_wrapper_no70beta175_errors, 'l013', 'var', /dobackground, /res36, /res25





makebg_gauss_iter, 4.25, 'l015', /doiter_var, /makelabel, /make36sizebg, /make25sizebg

gaussfit_iter_sedfit_wrapper_no70beta175_errors, 'l015', 'var', /dobackground, /res36, /res25





makebg_gauss_iter, 4.25, 'l017', /doiter_var, /makelabel, /make36sizebg, /make25sizebg

gaussfit_iter_sedfit_wrapper_no70beta175_errors, 'l017', 'var', /dobackground, /res36, /res25





makebg_gauss_iter, 4.25, 'l019', /doiter_var, /makelabel, /make36sizebg, /make25sizebg

gaussfit_iter_sedfit_wrapper_no70beta175_errors, 'l019', 'var', /dobackground, /res36, /res25




makebg_gauss_iter, 4.25, 'l021', /doiter_var, /makelabel, /make36sizebg, /make25sizebg

gaussfit_iter_sedfit_wrapper_no70beta175_errors, 'l021', 'var', /dobackground, /res36, /res25




makebg_gauss_iter, 4.25, 'l338', /doiter_var, /makelabel, /make36sizebg, /make25sizebg

gaussfit_iter_sedfit_wrapper_no70beta175_errors, 'l338', 'var', /dobackground, /res36, /res25




makebg_gauss_iter, 4.25, 'l341', /doiter_var, /makelabel, /make36sizebg, /make25sizebg

gaussfit_iter_sedfit_wrapper_no70beta175_errors, 'l341', 'var', /dobackground, /res36, /res25




makebg_gauss_iter, 4.25, 'l343', /doiter_var, /makelabel, /make36sizebg, /make25sizebg

gaussfit_iter_sedfit_wrapper_no70beta175_errors, 'l343', 'var', /dobackground, /res36, /res25




makebg_gauss_iter, 4.25, 'l345', /doiter_var, /makelabel, /make36sizebg, /make25sizebg

gaussfit_iter_sedfit_wrapper_no70beta175_errors, 'l345', 'var', /dobackground, /res36, /res25




makebg_gauss_iter, 4.25, 'l347', /doiter_var, /makelabel, /make36sizebg, /make25sizebg

gaussfit_iter_sedfit_wrapper_no70beta175_errors, 'l347', 'var', /dobackground, /res36, /res25




makebg_gauss_iter, 4.25, 'l349', /doiter_var, /makelabel, /make36sizebg, /make25sizebg

gaussfit_iter_sedfit_wrapper_no70beta175_errors, 'l349', 'var', /dobackground, /res36, /res25




makebg_gauss_iter, 4.25, 'l351', /doiter_var, /makelabel, /make36sizebg, /make25sizebg

gaussfit_iter_sedfit_wrapper_no70beta175_errors, 'l351', 'var', /dobackground, /res36, /res25




makebg_gauss_iter, 4.25, 'l354', /doiter_var, /makelabel, /make36sizebg, /make25sizebg

gaussfit_iter_sedfit_wrapper_no70beta175_errors, 'l354', 'var', /dobackground, /res36, /res25




makebg_gauss_iter, 4.25, 'l356', /doiter_var, /makelabel, /make36sizebg, /make25sizebg

gaussfit_iter_sedfit_wrapper_no70beta175_errors, 'l356', 'var', /dobackground, /res36, /res25




makebg_gauss_iter, 4.25, 'l358', /doiter_var, /makelabel, /make36sizebg, /make25sizebg

gaussfit_iter_sedfit_wrapper_no70beta175_errors, 'l358', 'var', /dobackground, /res36, /res25





EOF
