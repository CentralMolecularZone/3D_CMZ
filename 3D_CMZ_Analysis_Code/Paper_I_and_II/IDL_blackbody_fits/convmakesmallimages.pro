PRO convmakesmallimages, fwhm_want, newpixsize, prefix, $
                         regridimages=regridimages, convimages=convimages, $
                         inner40=inner40
                        

;ogfiles = array of 5 strings with data name (not directory) in this order
; sp500, sp350, sp250, pa170, pa70
; plw, pmw, psw, red, blue

; fwhm_want = 35.9
; newpixsize = 7.8
; prefix = 'l000'
; regridimages = 1
; convimages=1



;adjusted for NEW naming standard...  May 20, 2013
;ogfiles = [prefix+'/'+'destripe_'+prefix+'_PLW_wgls_rcal.fits', $
;           prefix+'/'+'destripe_'+prefix+'_PMW_wgls_rcal.fits', $
;           prefix+'/'+'destripe_'+prefix+'_PSW_wgls_rcal.fits', $
;           prefix+'/'+'destripe_'+prefix+'_red_wgls_rcal.fits', $
;           prefix+'/'+'destripe_'+prefix+'_blue_wgls_rcal.fits']


if keyword_set(inner40) then begin
;adjusted for inner40 naming standard...  March 6, 2013
ogfiles = [prefix+'/inner40_mosaic_500.fits', $
           prefix+'/inner40_mosaic_350.fits', $
           prefix+'/inner40_mosaic_250.fits', $
           prefix+'/inner40_mosaic_160.fits', $
           prefix+'/inner40_mosaic_70.fits']

endif else begin


ogfiles = strarr(5)
ogfiles(0) = file_search(prefix+'/HIGAL*500*')
ogfiles(1) = file_search(prefix+'/HIGAL*350*')
ogfiles(2) = file_search(prefix+'/HIGAL*250*')
ogfiles(3) = file_search(prefix+'/HIGAL*160*')
ogfiles(4) = file_search(prefix+'/HIGAL*70*')

endelse


; reffile = '/mira/scratch/cara/cloudcomp/l00/reffile_l00.fits'
if not keyword_set(fwhm_want) then fwhm_want = 35.9 & print, 'setting fwhm_want to 35.9'
if not keyword_set(newpixsize) then newpixsize = 7.8 & print, 'setting newpixsize to 7.8'

; fwhm_want = 35.9
; newpixsize = 7.8 
;    fwhm_want = 25.2
;    newpixsize= 5.5


fwhm_have = [35.9, 25.2, 18., 12.2, 5.]

if keyword_set(regridimages) then begin
    fits_read, ogfiles(0), sp500, sp500hdr
    fits_read, ogfiles(1), sp350, sp350hdr
    fits_read, ogfiles(2), sp250, sp250hdr
    fits_read, ogfiles(3), pa170, pa170hdr
    fits_read, ogfiles(4), pa70, pa70hdr

;HEULER and HASTROM all the files  
;fits_read, reffile, refdata, refhdr
    refhdr = pa70hdr            ;try this for now?
    
    heuler, sp500hdr, /galactic
    heuler, sp350hdr, /galactic
    heuler, sp250hdr, /galactic
    heuler, pa170hdr, /galactic
    heuler, pa70hdr, /galactic
    
    hastrom, sp500, sp500hdr, refhdr
    hastrom, sp350, sp350hdr, refhdr
    hastrom, sp250, sp250hdr, refhdr
    hastrom, pa170, pa170hdr, refhdr
    hastrom, pa70, pa70hdr, refhdr
    
    
;get rid of negative values before convolving!  set to NaN so they
;will be ignored in the convolution
    wneg = where( (sp500 le 0.) or (sp350 le 0.) or (sp250 le 0.) or $
                  (pa170 le 0.) or (pa70 le 0.) )
    sz = size(sp500)
    nanlabel = fltarr(sz(1), sz(2)) + 1
    nanlabel[wneg] = !values.f_nan
    
    sp500 *= nanlabel
    sp350 *= nanlabel
    sp250 *= nanlabel
    pa170 *= nanlabel
    pa70 *= nanlabel
    
    fits_write, prefix+'/'+'sp500_'+prefix+'_noconv.fits',$
      sp500, sp500hdr
    fits_write, prefix+'/'+'sp350_'+prefix+'_noconv.fits',$
      sp350, sp350hdr
    fits_write, prefix+'/'+'sp250_'+prefix+'_noconv.fits',$
      sp250, sp250hdr
    fits_write, prefix+'/'+'pa170_'+prefix+'_noconv.fits',$
      pa170,  pa170hdr
    fits_write, prefix+'/'+'pa70_'+prefix+'_noconv.fits',$
      pa70, pa70hdr
    
endif else begin
    
    fits_read, prefix+'/'+'sp500_'+prefix+'_noconv.fits',$
      sp500, sp500hdr
    fits_read, prefix+'/'+'sp350_'+prefix+'_noconv.fits',$
      sp350, sp350hdr
    fits_read, prefix+'/'+'sp250_'+prefix+'_noconv.fits',$
      sp250, sp250hdr
    fits_read, prefix+'/'+'pa170_'+prefix+'_noconv.fits',$
      pa170,  pa170hdr
    fits_read, prefix+'/'+'pa70_'+prefix+'_noconv.fits',$
      pa70, pa70hdr

endelse


if keyword_set(convimages) then begin

    sz = size(sp500)
    nanlabel = fltarr(sz(1), sz(2)) + 1
    wnan = where( finite(sp500, /nan) or finite(sp350, /nan) or $
                  finite(sp250, /nan) or $
                  finite(pa170, /nan) or finite(pa70, /nan) )
    nanlabel[wnan] = !values.f_nan

;Convolve to desired resolution using the new, improved convolve_res.pro!!!

    if fwhm_want gt 25.2 then sp350 = convolve_res(fwhm_want, fwhm_have(1), sp350, sp350hdr)
    
    sp250 = convolve_res(fwhm_want, fwhm_have(2), sp250, sp250hdr)
    pa170 = convolve_res(fwhm_want, fwhm_have(3), pa170, pa170hdr)
    pa70 = convolve_res(fwhm_want, fwhm_have(4), pa70, pa70hdr)

    sp500 *= nanlabel
    sp350 *= nanlabel
    sp250 *= nanlabel
    pa170 *= nanlabel
    pa70 *= nanlabel
    
;Regrid to common, reasonable pixel size
    pixx = 0 ; to make sure isn't saved from previous run
    pixx = abs(sxpar(pa70hdr, 'cd1_1'))*3600.
    if pixx eq 0 then pixx = abs(sxpar(pa70hdr, 'CDELT1'))*3600.
    print, 'pixx = ', pixx
   
    xx = newpixsize/pixx
    sz = size(sp500)
    szx = floor(sz(1)/xx)
    szy = floor(sz(2)/xx)
    
    hcongrid, sp500, sp500hdr, /half_half, outsize=[szx, szy]
    hcongrid, sp350, sp350hdr, /half_half, outsize=[szx, szy]
    hcongrid, sp250, sp250hdr, /half_half, outsize=[szx, szy]
    hcongrid, pa170, pa170hdr, /half_half, outsize=[szx, szy]
    hcongrid, pa70, pa70hdr, /half_half, outsize=[szx, szy]
    
    
;Write to new fits file
    fits_write, prefix+'/sp500_'+prefix+'_conv'+strn(round(fwhm_want))+'.fits', sp500, sp500hdr
    fits_write, prefix+'/sp350_'+prefix+'_conv'+strn(round(fwhm_want))+'.fits', sp350, sp350hdr
    fits_write, prefix+'/sp250_'+prefix+'_conv'+strn(round(fwhm_want))+'.fits', sp250, sp250hdr
    fits_write, prefix+'/pa170_'+prefix+'_conv'+strn(round(fwhm_want))+'.fits', pa170, pa170hdr
    fits_write, prefix+'/pa70_'+prefix+'_conv'+strn(round(fwhm_want))+'.fits', pa70, pa70hdr


endif







END

