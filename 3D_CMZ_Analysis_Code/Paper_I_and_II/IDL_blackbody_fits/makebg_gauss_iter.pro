;PRO makebg_l00_gauss_iter, limitnum, cv, doiter=doiter, makelabel=makelabel, $
;                           make36sizebg=make36sizebg, make25sizebg=make25sizebg

 PRO makebg_gauss_iter, limitnum, prefix, doiter_fixed=doiter_fixed, $
                        doiter_orig=doiter_orig, doiter_var=doiter_var, $
                        makelabel=makelabel, $
                        make36sizebg=make36sizebg, make25sizebg=make25sizebg
 
 
;limitnum is either the X*sigma for doiter_orig and doiter_var
; and is the definite flux limit for doiter_fixed                    

;doiter_orig does the iterations as outlined in Battersby et al. 2011
; does 20 iterations and chooses iteration 16

;doiter_var does iterations until ...

;doiter_fixed sets a specific cutoff 500 micron flux value and creates
; the label map etc. for that cutoff
; set to the value of the 500 micron cutoff!

;;;
;;;run with cv = '36' only !!!
;;;

; limitnum = 4.25
; prefix = 'l033'
; doiter_var = 1
; makelabel = 1
; make36sizebg = 1
; make25sizebg = 1

; cv = '36'
; doiter=1
; makelabel=1
; make25sizebg=1
; make36sizebg=1

print, 'limitnum', limitnum

limstr=strn(limitnum)

;;;
;;;run with cv = '36' only !!!
;;;

cv = '36'
cvfloat = float(cv)
cvstr = cv

fits_read, prefix+'/sp500_'+prefix+'_conv'+cv+'.fits', sp500, sp500hdr
fits_read, prefix+'/sp350_'+prefix+'_conv'+cv+'.fits', sp350, sp350hdr
fits_read, prefix+'/sp250_'+prefix+'_conv'+cv+'.fits', sp250, sp250hdr
fits_read, prefix+'/pa170_'+prefix+'_conv'+cv+'.fits', pa170, pa170hdr
fits_read, prefix+'/pa70_'+prefix+'_conv'+cv+'.fits', pa70, pa70hdr

   
;Get coordinates and save to b.sav for gauss fit across plane
centerx = sxpar(sp500hdr, 'CRPIX1')
centery = sxpar(sp500hdr, 'CRPIX2')
naxis1 = sxpar(sp500hdr, 'NAXIS1') 
naxis2 = sxpar(sp500hdr, 'NAXIS2')
centerl = sxpar(sp500hdr, 'CRVAL1')
centerb = sxpar(sp500hdr, 'CRVAL2')

centerx = centerx-1      ;idl, ds9
centery = centery-1

x0 = centerx - (naxis1/2)+1
x1 = centerx + (naxis1/2)
y0 = centery - (naxis2/2)+1
y1 = centery + (naxis2/2)

dx = 0 & dy = 0
dx = sxpar(sp500hdr, 'CDELT1') & dy = sxpar(sp500hdr, 'CDELT2')
if dx eq 0 then begin
dx = sxpar(sp500hdr, 'cd1_1')
dy = sxpar(sp500hdr, 'cd2_2')
endif
if (dx eq 0) or (dx eq 1) then begin
dx = sxpar(sp500hdr, 'PC1_1')
dy = sxpar(sp500hdr, 'PC2_2')
endif

if dx eq 0 then dx = sxpar(sp500hdr, 'CDELT1')
if dy eq 0 then dy = sxpar(sp500hdr, 'CDELT2')

centerx = naxis1/2
centery = naxis2/2

x = lindgen(naxis1)
y = lindgen(naxis2)
l = (x-centerx)*dx + centerl
if max(l) gt 360 then l-=360
b = (y-centery)*dy+centerb

save, b, filename=prefix+'/b_'+prefix+'.sav'


sz = size(sp350)

wnan = where( finite(sp500, /nan) or finite(sp350, /nan) or finite(sp250, /nan) or $
              finite(pa170, /nan) )


if keyword_set(doiter_orig) then begin
print, 'doiter_orig set!  Proceeding ...'
;Iterate on gaussian fit to background
    niter = 20 ;on _small file, iterations converged around iteration 8
;to create the first bg image for comparison
    firstbg = 1
    origfile = prefix+'/sp500_'+prefix+'_conv'+cv+'.fits'
    fits_read, origfile, ogdata, oghdr
    if keyword_set(firstbg) then begin
        bg0 = convolve_res(720., cvfloat, ogdata, oghdr)
        bg0[wnan] = !values.f_nan
        fits_write, prefix+'/bgimage0_conv720_'+prefix+'_gaussfit.fits', bg0, oghdr
    endif 
    
    limitarr = fltarr(niter)

    for ii = 0, niter-1 do begin
        iterstr = strn(ii+1) ;convert to string and remove padded blanks
        iterstrmin1 = strn(ii) ;bg image used is from the previous iteration
        origfile = prefix+'/sp500_'+prefix+'_conv'+cvstr+'.fits'
        bgfile = prefix+'/bgimage'+iterstrmin1+'_conv720_'+prefix+'_gaussfit.fits'
        gauss_iter_wrapper, origfile, bgfile, $
          prefix, b, cvstr, iterstr, limitnum, wnan, limit
        limitarr[ii] = limit
    endfor
limits = limitarr
save, limits, limitarr, filename=prefix+'/limit_'+prefix+'_orig_'+limstr+'sigma.sav'    
endif

;decided almost arbitrarily, but that's okay!
iterstr = '16'
iterstrmin1 = '15'


if keyword_set(doiter_var) then begin
print, 'doiter_var set!  Proceeding ...'
;Iterate on gaussian fit to background
                                ; do at least 10 iterations
    niter = 3.
    firstbg = 1
    origfile = prefix+'/sp500_'+prefix+'_conv'+cv+'.fits'
    fits_read, origfile, ogdata, oghdr
    if keyword_set(firstbg) then begin
        bg0 = convolve_res(720., cvfloat, ogdata, oghdr)
        bg0[wnan] = !values.f_nan
        fits_write, prefix+'/bgimage0_conv720_'+prefix+'_gaussfit.fits', bg0, oghdr
    endif 
    
    limitarr = fltarr(100)      ;give plenty of room to start
    
                                ;do at least niter (10) iterations to start
    for ii = 0, niter-1 do begin
        iterstr = strn(ii+1) ;convert to string and remove padded blanks
        iterstrmin1 = strn(ii) ;bg image used is from the previous iteration
        origfile = prefix+'/sp500_'+prefix+'_conv'+cvstr+'.fits'
        bgfile = prefix+'/bgimage'+iterstrmin1+'_conv720_'+prefix+'_gaussfit.fits'
        gauss_iter_wrapper, origfile, bgfile, $
          prefix, b, cvstr, iterstr, limitnum, wnan, limit
        limitarr[ii] = limit
    endfor
    
;now keep iterating until ... this plus the last 2 values has a standard
;deviation of less than 1
    sval = 100.
    
    while sval gt 1 do begin
        ii = ii +1
        iterstr = strn(ii) ;convert to string and remove padded blanks
        iterstrmin1 = strn(ii-1) ;bg image used is from the previous iteration
        origfile = prefix+'/sp500_'+prefix+'_conv'+cvstr+'.fits'
        bgfile = prefix+'/bgimage'+iterstrmin1+'_conv720_'+prefix+'_gaussfit.fits'
        gauss_iter_wrapper, origfile, bgfile, $
          prefix, b, cvstr, iterstr, limitnum, wnan, limit
        limitarr[ii] = limit
        
        lmvals = [limitarr(ii), limitarr(ii-1), limitarr(ii-2)]
        sval = stdev(lmvals)

        ;ii = ii+1
        
    endwhile
    
   
    limits = limitarr

    save, limits, limitarr, filename=prefix+'/limit_'+prefix+'_var_'+limstr+'sigma.sav'    


    
;decided NOT arbitrarily, so that's great!
    iterstr = strn(ii)
    iterstrmin1 = strn(ii-1)
    print, 'final iteration: ', iterstr
    print, 'final limit: ', limitarr(ii)

endif



if keyword_set(doiter_fixed) then begin
print, 'doiter_fixed set!  Proceeding ...'
;Iterate on gaussian fit to background
;    niter = 20 ;on _small file, iterations converged around iteration 8
;to create the first bg image for comparison
    firstbg = 1
    origfile = prefix+'/sp500_'+prefix+'_conv'+cv+'.fits'
    fits_read, origfile, ogdata, oghdr
    if keyword_set(firstbg) then begin
        bg0 = convolve_res(720., cvfloat, ogdata, oghdr)
        bg0[wnan] = !values.f_nan
        fits_write, prefix+'/bgimage0_conv720_'+prefix+'_gaussfit.fits', bg0, oghdr
    endif 
    
    origfile = prefix+'/sp500_'+prefix+'_conv'+cvstr+'.fits'
    bgfile = prefix+'/bgimage0_conv720_'+prefix+'_gaussfit.fits' 
    gauss_iter_wrapper_fixedval, origfile, bgfile, prefix, b, cv, wnan, limitnum
    iterstr='fixed'
    limits = limitnum

    save, limits, filename=prefix+'/limit_'+prefix+'_fixed_'+limstr+'.sav'    


endif



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Make Backgrounds at Every Wavelength ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if keyword_set(makelabel) then begin
    print, 'making label'
    fits_read, prefix+'/bgimage'+iterstr+'_conv720_'+prefix+'_gaussfit.fits', bg500, bg500hdr
    fits_read, prefix+'/gfit_label_iter'+iterstr+'_'+prefix+'.fits', label8, label8hdr
    
;get rid of small regions in label map (less than 10 pixels or so)
    sz = size(label8)
    label = label8
    
    for ii = 2, sz(1) - 3 do begin
        for jj = 2, sz(2) - 3 do begin
                                ;make 5x5 box around pixel
            area = 5.*5.
            test = label8[ii-2:ii+2, jj-2:jj+2]
            totalmask = total(test)
            if totalmask lt 10 then label[ii,jj] = 0 else label[ii,jj] = label8[ii,jj]
        endfor
    endfor    

    qq = where(label eq 1, nmask)
;    finallimit = limits(15)
;    openw, 1, 'sigmatest.txt', /append
;    printf, 1, limitnum, nmask, finallimit
;    close, 1

if keyword_set(doiter_orig) then iterstr='orig'
if keyword_set(doiter_var) then iterstr='var'
if keyword_set(doiter_fixed) then iterstr='fixed'

    fits_write, prefix+'/gfit_iterations_label_'+prefix+'_iter'+iterstr+'.fits', label, oghdr
    
endif else fits_read, prefix+'/gfit_iterations_label_'+prefix+'_iter'+iterstr+'.fits', label, labelhdr



if keyword_set(make36sizebg) then begin
nanlabel = label
wm = where(label eq 1, complement=wn) ;wm = where mask, wn = where not mask
;masks are nans, else, 1's.  this is now an image where the masks are
;nan's and everything else is one... so it can be multiplied by any
;image to create a "masked" image
nanlabel[wm] = !values.f_nan  
nanlabel[wn] = 1



print, 'making backgrounds'

;SPIRE 500 microns
sp500masked = sp500*nanlabel
sp500bg = convolve_res(720., cvfloat, sp500masked, sp500hdr)
sp500bg[wnan] = !values.f_nan
fits_write, prefix+'/sp500bg_'+prefix+'_iter'+iterstr+'_gaussfit.fits', sp500bg, sp500hdr
;inside the mask, it is the difference of the original image and the
;background, outside, simply the background
sp500diff = sp500
sp500diff[wm] = sp500[wm] - sp500bg[wm]
sp500diff[wn] = sp500bg[wn]
fits_write, prefix+'/sp500diff_'+prefix+'_iter'+iterstr+'_gaussfit.fits', sp500diff, sp500hdr

;SPIRE 350 microns
sp350masked = sp350*nanlabel
sp350bg = convolve_res(720., cvfloat, sp350masked, sp350hdr)
sp350bg[wnan] = !values.f_nan
fits_write, prefix+'/sp350bg_'+prefix+'_iter'+iterstr+'_gaussfit.fits', sp350bg, sp500hdr
;inside the mask, it is the difference of the original image and the
;background, outside, simply the background
sp350diff = sp350
sp350diff[wm] = sp350[wm] - sp350bg[wm]
sp350diff[wn] = sp350bg[wn]
fits_write, prefix+'/sp350diff_'+prefix+'_iter'+iterstr+'_gaussfit.fits', sp350diff, sp500hdr

;SPIRE 250 microns
sp250masked = sp250*nanlabel
sp250bg = convolve_res(720., cvfloat, sp250masked, sp250hdr)
sp250bg[wnan] = !values.f_nan
fits_write, prefix+'/sp250bg_'+prefix+'_iter'+iterstr+'_gaussfit.fits', sp250bg, sp500hdr
;inside the mask, it is the difference of the original image and the
;background, outside, simply the background
sp250diff = sp250
sp250diff[wm] = sp250[wm] - sp250bg[wm]
sp250diff[wn] = sp250bg[wn]
fits_write, prefix+'/sp250diff_'+prefix+'_iter'+iterstr+'_gaussfit.fits', sp250diff, sp500hdr

;PACS 170 microns
pa170masked = pa170*nanlabel
pa170bg = convolve_res(720., cvfloat, pa170masked, pa170hdr)
pa170bg[wnan] = !values.f_nan
fits_write, prefix+'/pa170bg_'+prefix+'_iter'+iterstr+'_gaussfit.fits', pa170bg, sp500hdr
;inside the mask, it is the difference of the original image and the
;background, outside, simply the background
pa170diff = pa170
pa170diff[wm] = pa170[wm] - pa170bg[wm]
pa170diff[wn] = pa170bg[wn]
fits_write, prefix+'/pa170diff_'+prefix+'_iter'+iterstr+'_gaussfit.fits', pa170diff, sp500hdr

;PACS 70 microns
pa70masked = pa70*nanlabel
pa70bg = convolve_res(720., cvfloat, pa70masked, pa70hdr)
pa70bg[wnan] = !values.f_nan
fits_write, prefix+'/pa70bg_'+prefix+'_iter'+iterstr+'_gaussfit.fits', pa70bg, sp500hdr
;inside the mask, it is the difference of the original image and the
;background, outside, simply the background
pa70diff = pa70
pa70diff[wm] = pa70[wm] - pa70bg[wm]
pa70diff[wn] = pa70bg[wn]
fits_write, prefix+'/pa70diff_'+prefix+'_iter'+iterstr+'_gaussfit.fits', pa70diff, sp500hdr

endif


if keyword_set(make25sizebg) then begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Make Backgrounds at 25" resolution ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

fits_read, prefix+'/sp500_'+prefix+'_conv25.fits', $
   sp500_25, sp500_25hdr
fits_read, prefix+'/sp350_'+prefix+'_conv25.fits', $
   sp350_25, sp350_25hdr
fits_read, prefix+'/sp250_'+prefix+'_conv25.fits', $
   sp250_25, sp250_25hdr
fits_read, prefix+'/pa170_'+prefix+'_conv25.fits', $
   pa170_25, pa170_25hdr
fits_read, prefix+'/pa70_'+prefix+'_conv25.fits', $
   pa70_25, pa70_25hdr

sz = size(sp500_25)

wnan = where( finite(sp500_25, /nan) or finite(sp350_25, /nan) or finite(sp250_25, /nan) or $
              finite(pa170_25, /nan) )

;regrid label plot to 25" pixel grid
fits_read, prefix+'/gfit_iterations_label_'+prefix+'_iter'+iterstr+'.fits', label36, label36hdr
hcongrid, label36, label36hdr, label25, label25hdr, /half, outsize=[sz(1), sz(2)]
fits_write, prefix+'/gfit_iterations_label_'+prefix+'_iter'+iterstr+'_conv25.fits', label25, label25hdr

 nanlabel25 = label25
 wm25 = where(label25 gt 0, complement=wn25) ;wm = where mask, wn = where not mask
 ;masks are nans, else, 1's.  this is now an image where the masks are
 ;nan's and everything else is one... so it can be multiplied by any
 ;image to create a "masked" image
 nanlabel25[wm25] = !values.f_nan  
 nanlabel25[wn25] = 1


;SPIRE 500 microns
;regrid bg to 25" pixel grid
fits_read, prefix+'/sp500bg_'+prefix+'_iter'+iterstr+'_gaussfit.fits', sp500bg36, sp500bg36hdr
hcongrid, sp500bg36, sp500bg36hdr, sp500bg25, sp500bg25hdr, /half, outsize=[sz(1), sz(2)]
sp500diff25 = sp500_25
sp500diff25[wm25] = sp500_25[wm25] - sp500bg25[wm25]
sp500diff25[wn25] = sp500bg25[wn25]
sp500diff25[wnan] = !values.f_nan
fits_write, prefix+'/sp500diff_'+prefix+'_iter'+iterstr+'_gaussfit_conv25.fits', sp500diff25, sp500_25hdr

;SPIRE 350 microns
;regrid bg to 25" pixel grid
fits_read, prefix+'/sp350bg_'+prefix+'_iter'+iterstr+'_gaussfit.fits', sp350bg36, sp350bg36hdr
hcongrid, sp350bg36, sp350bg36hdr, sp350bg25, sp350bg25hdr, /half, outsize=[sz(1), sz(2)]
sp350diff25 = sp350_25
sp350diff25[wm25] = sp350_25[wm25] - sp350bg25[wm25]
sp350diff25[wn25] = sp350bg25[wn25]
sp350diff25[wnan] = !values.f_nan
fits_write, prefix+'/sp350diff_'+prefix+'_iter'+iterstr+'_gaussfit_conv25.fits', sp350diff25, sp350_25hdr

;SPIRE 250 microns
;regrid bg to 25" pixel grid
fits_read, prefix+'/sp250bg_'+prefix+'_iter'+iterstr+'_gaussfit.fits', sp250bg36, sp250bg36hdr
hcongrid, sp250bg36, sp250bg36hdr, sp250bg25, sp250bg25hdr, /half, outsize=[sz(1), sz(2)]
sp250diff25 = sp250_25
sp250diff25[wm25] = sp250_25[wm25] - sp250bg25[wm25]
sp250diff25[wn25] = sp250bg25[wn25]
sp250diff25[wnan] = !values.f_nan
fits_write, prefix+'/sp250diff_'+prefix+'_iter'+iterstr+'_gaussfit_conv25.fits', sp250diff25, sp250_25hdr

;PACS 170 microns
;regrid bg to 25" pixel grid
fits_read, prefix+'/pa170bg_'+prefix+'_iter'+iterstr+'_gaussfit.fits', pa170bg36, pa170bg36hdr
hcongrid, pa170bg36, pa170bg36hdr, pa170bg25, pa170bg25hdr, /half, outsize=[sz(1), sz(2)]
pa170diff25 = pa170_25
pa170diff25[wm25] = pa170_25[wm25] - pa170bg25[wm25]
pa170diff25[wn25] = pa170bg25[wn25]
pa170diff25[wnan] = !values.f_nan
fits_write, prefix+'/pa170diff_'+prefix+'_iter'+iterstr+'_gaussfit_conv25.fits', pa170diff25, pa170_25hdr

;PACS 70 microns
;regrid bg to 25" pixel grid
fits_read, prefix+'/pa70bg_'+prefix+'_iter'+iterstr+'_gaussfit.fits', pa70bg36, pa70bg36hdr
hcongrid, pa70bg36, pa70bg36hdr, pa70bg25, pa70bg25hdr, /half, outsize=[sz(1), sz(2)]
pa70diff25 = pa70_25
pa70diff25[wm25] = pa70_25[wm25] - pa70bg25[wm25]
pa70diff25[wn25] = pa70bg25[wn25]
pa70diff25[wnan] = !values.f_nan
fits_write, prefix+'/pa70diff_'+prefix+'_iter'+iterstr+'_gaussfit_conv25.fits', pa70diff25, pa70_25hdr


endif



END
