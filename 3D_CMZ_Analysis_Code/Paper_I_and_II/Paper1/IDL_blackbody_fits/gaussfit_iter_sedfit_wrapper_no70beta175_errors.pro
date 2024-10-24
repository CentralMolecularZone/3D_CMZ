
;No longer has option to include BGPS, BGPS*1.5 or to include the 70
;micron point in the dense clump fitting.
; The method now is to include SPIRE 500, 350, 250, and PACS 170 and
; 70 microns in the background fitting and all the above except PACS
; 70 in the dense clump fitting.
;beta is set to 1.75 in the dense clump fitting, and is left as a free
;parameter in the background fitting.
PRO gaussfit_iter_sedfit_wrapper_no70beta175_errors, prefix, iterstr, $
                                                     dobackground=dobackground, $
                                                     res36=res36, res25=res25
;iterstr options
;iterstr = 'var'
;iterstr = 'orig'
;iterstr = 'fixed'


; dobackground = 0
; res36 = 0
; res25 = 1
; prefix = 'l000'
; dobgps = 0
; dobgps15 = 0
; no70 = 1
;beta175 = 1
;beta2 = 0

; if keyword_set(res36) then cv='36'
; if keyword_set(res25) then cv='25'

; if not keyword_set(cv) then begin
;     print, '!!! need to set res36 or res25 !!!'
;     stop
; endif

;halfstr = 'highhalf'

; to start for bg stuff ...
cv = '36'

;not from 4sigma directory

;from makebg_gauss_iter.pro
fits_read, prefix+'/gfit_iterations_label_'+prefix+'_iter'+iterstr+'.fits', label, labelhdr
fits_read, prefix+'/gfit_iterations_label_'+prefix+'_iter'+iterstr+'_conv25.fits', label25, label25hdr

print, 'here1'
;Raw, convolved, regridded, MJy/Sr files
fits_read, prefix+'/sp500_'+prefix+'_conv'+cv+'.fits', $
  sp500, sp500hdr
fits_read, prefix+'/sp350_'+prefix+'_conv'+cv+'.fits', $
  sp350, sp350hdr
fits_read, prefix+'/sp250_'+prefix+'_conv'+cv+'.fits', $
  sp250, sp250hdr
fits_read, prefix+'/pa170_'+prefix+'_conv'+cv+'.fits', $
  pa170, pa170hdr
fits_read, prefix+'/pa70_'+prefix+'_conv'+cv+'.fits', $
  pa70, pa70hdr

print, 'here2'
;Smoothed background images (from makebg_gauss_iter.pro)
fits_read, prefix+'/sp500bg_'+prefix+'_iter'+iterstr+'_gaussfit.fits', bg500, bg500hdr
fits_read, prefix+'/sp350bg_'+prefix+'_iter'+iterstr+'_gaussfit.fits', bg350, bg350hdr
fits_read, prefix+'/sp250bg_'+prefix+'_iter'+iterstr+'_gaussfit.fits', bg250, bg250hdr
fits_read, prefix+'/pa170bg_'+prefix+'_iter'+iterstr+'_gaussfit.fits', bg170, bg170hdr
fits_read, prefix+'/pa70bg_'+prefix+'_iter'+iterstr+'_gaussfit.fits', bg70, bg70hdr

print, 'here3'
;Data - BG, masked to bolocat points (outside that = 0)
; from makebg_gauss_iter.pro)
;CONV36
fits_read, prefix+'/sp500diff_'+prefix+'_iter'+iterstr+'_gaussfit.fits',$
  sp500diffmask, sp500diffmaskhdr
fits_read, prefix+'/sp350diff_'+prefix+'_iter'+iterstr+'_gaussfit.fits',$
  sp350diffmask, sp350diffmaskhdr
fits_read, prefix+'/sp250diff_'+prefix+'_iter'+iterstr+'_gaussfit.fits',$
  sp250diffmask, sp250diffmaskhdr
fits_read, prefix+'/pa170diff_'+prefix+'_iter'+iterstr+'_gaussfit.fits',$
  pa170diffmask, pa170diffmaskhdr
fits_read, prefix+'/pa70diff_'+prefix+'_iter'+iterstr+'_gaussfit.fits',$
  pa70diffmask, pa70diffmaskhdr

wz = where(label eq 0, complement=ws)

print, 'here4'
c = 2.99792458e8 ;m/s

npix = n_elements(sp350)
sz = size(sp350)
temp = fltarr(sz(1), sz(2))
column = fltarr(sz(1), sz(2))
beta = fltarr(sz(1), sz(2))
tempbg = fltarr(sz(1), sz(2))
columnbg = fltarr(sz(1), sz(2))
betabg = fltarr(sz(1), sz(2))
label36 = fltarr(sz(1), sz(2))
bgchisqr = fltarr(sz(1), sz(2))
bgdegfree = fltarr(sz(1), sz(2))
bgtemperr = fltarr(sz(1), sz(2))
bgcolerr = fltarr(sz(1), sz(2))
bgbetaerr = fltarr(sz(1), sz(2))

;FREE BETA
  parinfo_free = replicate({value:0.D, limited:[0,0], limits:[0.D,0]},3)

  parinfo_free[2].limited[0:1] = [1,1]
  parinfo_free[2].limits[0:1] = [1,3]
  parinfo_free[*].value = [20., 1., 1.75]

  parinfo_free[0].limited[0:1] = [1,1]
  parinfo_free[0].limits[0:1] = [0,100] ;temp

  parinfo_free[1].limited[0:1] = [1,1]
  parinfo_Free[1].limits[0:1] = [0,1000]  ;column, max now at 1000


;FIXED BETA
 parinfo_fix2 = replicate({value:0.D, limited:[0,0], limits:[0.D,0]},2)
 parinfo_fix2[*].value = [20., 1.]

 parinfo_fix2[0].limited[0:1] = [1,1]
 parinfo_fix2[0].limits[0:1] = [0,100] ;temp

 parinfo_fix2[1].limited[0:1] = [1,1]
 parinfo_fix2[1].limits[0:1] = [0,1000]  ;column, max now at 1000


;lambda = [1.1101e-3, 520e-6, 360e-6, 250e-6, 160e-6, 70e-6]
lambda = [520e-6, 360e-6, 250e-6, 160e-6, 70e-6]
nu = c/lambda
nuall = fill_array(1000, min(nu), max(nu))
np = n_elements(nu) -1 ;number of points-1 = index of last value

;for ignoring nan points when necessary


;BACKGROUND 
if dobackground then begin
nubg = c/lambda
    for jj = 0, sz(2)-1 do begin
        for ii = 0, sz(1)-1 do begin
;   ii = 121
;   jj = 903
; print, 'here'

;OLD STRATEGY, IGNORE THOSE POINTS
;             if (bg500(ii,jj) le 0.) or $
;               (bg350(ii,jj) le 0.) or (bg250(ii,jj) le 0.) or $
;               (bg170(ii,jj) le 0.) or (bg70(ii,jj) le 0.) then begin 
;                 tempbg[ii,jj] = 0.
;                 columnbg[ii,jj] = 0.
;                 betabg[ii,jj] = 0.
;             endif else begin

;got rid of negatives.
;because only have negative values remaining at PACS 170....
; if bg170(ii,jj) le 60. then begin
;if finite(nan170(ii,jj), /nan) eq 1 then begin
if bg170(ii,jj) le bg250(ii,jj) then begin
     inubg = [bg500(ii,jj), bg350(ii,jj), bg250(ii,jj), bg70(ii,jj)]
     err = inubg*0.2
     lambdabg = [520e-6, 360e-6, 250e-6,  70e-6]
     nubg = c/lambdabg
 endif else begin



                inubg = [bg500(ii,jj),bg350(ii,jj),$
                         bg250(ii,jj),bg170(ii,jj),bg70(ii,jj)]
                err = inubg*0.2
                nubg = c/lambda
endelse  

if finite(bg70(ii,jj), /nan) or finite(bg170(ii,jj), /nan) or $
  finite(bg250(ii,jj), /nan) or finite(bg350(ii,jj), /nan) or $
  finite(bg500(ii,jj), /nan) then begin
    
    tempbg[ii,jj] = !values.f_nan
    columnbg[ii,jj] = !values.f_nan
    betabg[ii,jj] = !values.f_nan
    bgtemperr[ii,jj] = !values.f_nan
    bgcolerr[ii,jj] = !values.f_nan
    bgbetaerr[ii,jj] = !values.f_nan

endif else begin

                pfitbg = mpfitfun('greybody_freebeta', nubg, inubg, err, $
                                  /quiet, parinfo=parinfo_free, bestnorm=chi, dof=freedeg, covar=comatrix)
                
                bgchisqr[ii,jj] = chi
                bgdegfree[ii,jj] = freedeg
                bgtemperr[ii,jj] = sqrt(comatrix(0,0))
                bgcolerr[ii,jj] = sqrt(comatrix(1,1))
                bgbetaerr[ii,jj] = sqrt(comatrix(2,2))
                
                tempbg[ii,jj] = pfitbg(0)
                columnbg[ii,jj] = pfitbg(1)
                betabg[ii,jj] = pfitbg(2)
endelse


        endfor
    endfor

wnan = where( finite(sp500, /nan) or finite(sp350, /nan) or finite(sp250, /nan) or $
              finite(pa170, /nan) )

tempbg[wnan] = !values.f_nan
columnbg[wnan]= !values.f_nan
betabg[wnan] = !values.f_nan
bgtemperr[wnan] = !values.f_nan
bgcolerr[wnan] = !values.f_nan
bgchisqr[wnan] = !values.f_nan
bgbetaerr[wnan] = !values.f_nan

fits_write, prefix+'/bg_temp_'+prefix+'_iter'+iterstr+'_gaussfit_iter_conv36.fits', tempbg, sp500hdr
fits_write, prefix+'/bg_column_'+prefix+'_iter'+iterstr+'_gaussfit_iter_conv36.fits', columnbg, sp500hdr
fits_write, prefix+'/bg_beta_'+prefix+'_iter'+iterstr+'_gaussfit_iter_conv36.fits', betabg, sp500hdr
fits_write, prefix+'/bg_temperr_'+prefix+'_iter'+iterstr+'_gaussfit_iter_conv36.fits', bgtemperr, sp500hdr
fits_write, prefix+'/bg_colerr_'+prefix+'_iter'+iterstr+'_gaussfit_iter_conv36.fits', bgcolerr, sp500hdr
fits_write, prefix+'/bg_chisqr_'+prefix+'_iter'+iterstr+'_gaussfit_iter_conv36.fits', bgchisqr, sp500hdr
fits_write, prefix+'/bg_betaerr_'+prefix+'_iter'+iterstr+'_gaussfit_iter_conv36.fits', bgbetaerr, sp500hdr

endif else begin
 
fits_read, prefix+'/bg_temp_'+prefix+'_iter'+iterstr+'_gaussfit_iter_conv36.fits', tempbg, hdr
fits_read, prefix+'/bg_column_'+prefix+'_iter'+iterstr+'_gaussfit_iter_conv36.fits', columnbg, hdr
fits_read, prefix+'/bg_beta_'+prefix+'_iter'+iterstr+'_gaussfit_iter_conv36.fits', betabg, hdr
fits_read, prefix+'/bg_temperr_'+prefix+'_iter'+iterstr+'_gaussfit_iter_conv36.fits', bgtemperr, hdr
fits_read, prefix+'/bg_colerr_'+prefix+'_iter'+iterstr+'_gaussfit_iter_conv36.fits', bgcolerr, hdr
fits_read, prefix+'/bg_chisqr_'+prefix+'_iter'+iterstr+'_gaussfit_iter_conv36.fits', bgchisqr, hdr
fits_read, prefix+'/bg_betaerr_'+prefix+'_iter'+iterstr+'_gaussfit_iter_conv36.fits', bgbetaerr, hdr

endelse


chisqr36 = fltarr(sz(1), sz(2))
degfree36 = fltarr(sz(1), sz(2))
temperr36 = fltarr(sz(1), sz(2))
colerr36 = fltarr(sz(1), sz(2))

nb = 0
nn = np-1

;Main image (36" resolution)
if keyword_set(res36) then begin
cv='36'
for jj = 0, sz(2)-1 do begin
    for ii = 0, sz(1)-1 do begin
        if label[ii,jj] eq 0 then begin
            temp[ii,jj] = tempbg(ii,jj)
            column[ii,jj] = columnbg(ii,jj)
            label36[ii,jj] = 0
            chisqr36[ii,jj] = bgchisqr(ii,jj)
            temperr36[ii,jj] = bgtemperr(ii,jj)
            colerr36[ii,jj] = bgcolerr(ii,jj)
        endif else begin
            if (sp500diffmask(ii,jj) le 0.) or $
              (sp350diffmask(ii,jj) le 0.) or (sp250diffmask(ii,jj) le 0.) or $
              (pa170diffmask(ii,jj) le 0.) or finite(pa170diffmask(ii,jj), /nan) or $
              finite(sp250diffmask(ii,jj), /nan) or finite(sp350diffmask(ii,jj), /nan) or $
              finite(sp500diffmask(ii,jj), /nan) then begin
                
                temp[ii,jj] = tempbg(ii,jj)
                column[ii,jj] = columnbg(ii,jj)
                label36[ii,jj] = 0
                chisqr36[ii,jj] = bgchisqr(ii,jj)
                temperr36[ii,jj] = bgtemperr(ii,jj)
                colerr36[ii,jj] = bgcolerr(ii,jj)
            endif else begin
                inu = [sp500diffmask(ii,jj),sp350diffmask(ii,jj),$
                       sp250diffmask(ii,jj),pa170diffmask(ii,jj),pa70diffmask(ii,jj)]
                err = inu*0.2
                inu = inu[nb:nn]
                err = err[nb:nn]
                nun = nu[nb:nn]
                pfit = mpfitfun('greybody_beta175', nun, inu, err, $
                                /quiet, parinfo=parinfo_fix2, bestnorm=chi, dof=freedeg, covar=comatrix)
    
                chisqr36[ii,jj] = chi
                degfree36[ii,jj] = freedeg
                temperr36[ii,jj] = sqrt(comatrix(0,0))
                colerr36[ii,jj] = sqrt(comatrix(1,1))

                temp[ii,jj] = pfit(0)
                column[ii,jj] = pfit(1)
                label36[ii,jj] = 1
            endelse
        endelse
;  plot, nun, inu, /xlog, /ylog, psym=1
;  oplot, nuall, greybody_beta175(nuall, pfit)
; print, pfit

; stop
    endfor
endfor


wnan = where( finite(sp500, /nan) or finite(sp350, /nan) or finite(sp250, /nan) or $
              finite(pa170, /nan) )

temp[wnan] = !values.f_nan
column[wnan]= !values.f_nan
label36[wnan] = !values.f_nan
temperr36[wnan] = !values.f_nan
colerr36[wnan] = !values.f_nan
chisqr36[wnan] = !values.f_nan

fits_write, prefix+'/gaussfit_iter_beta175_temp_'+prefix+'_iter'+iterstr+'_conv36.fits', temp, sp500hdr
fits_write, prefix+'/gaussfit_iter_beta175_column_'+prefix+'_iter'+iterstr+'_conv36.fits', column, sp500hdr
fits_write, prefix+'/gaussfit_iter_beta175_column_'+prefix+'_iter'+iterstr+'_conv36_label.fits', label36, sp500hdr
fits_write, prefix+'/gaussfit_iter_beta175_temperr_'+prefix+'_iter'+iterstr+'_conv36.fits', temperr36, sp500hdr
fits_write, prefix+'/gaussfit_iter_beta175_colerr_'+prefix+'_iter'+iterstr+'_conv36.fits', colerr36, sp500hdr
fits_write, prefix+'/gaussfit_iter_beta175_chisqr_'+prefix+'_iter'+iterstr+'_conv36.fits', chisqr36, sp500hdr


endif


;;;;;;;;;;;;;;;;;;;;;;
;;; 25" resolution ;;;
;;;;;;;;;;;;;;;;;;;;;;
if keyword_set(res25) then begin
cv='25'

nb = 1 ;no sp500 point

;LABELMAP file with ones
fits_read, prefix+'/gfit_iterations_label_'+prefix+'_iter'+iterstr+'_conv25.fits', $
  label, labelhdr

;Raw, convolved, regridded, MJy/Sr files
; from convmakesmallimages_'+prefix+'highhalf.pro and convmakesmallimages_'+prefix+'lowhalf.pro
fits_read, prefix+'/sp500_'+prefix+''+'_conv25.fits', $
  sp500, sp500hdr
fits_read, prefix+'/sp350_'+prefix+''+'_conv25.fits', $
  sp350, sp350hdr
fits_read, prefix+'/sp250_'+prefix+''+'_conv25.fits', $
  sp250, sp250hdr
fits_read, prefix+'/pa170_'+prefix+''+'_conv25.fits', $
  pa170, pa170hdr
fits_read, prefix+'/pa70_'+prefix+''+'_conv25.fits', $
  pa70, pa70hdr


;Smoothed background images (from makebg_'+prefix+'_gauss_iter.pro)
fits_read, prefix+'/sp500bg_'+prefix+'_iter'+iterstr+'_gaussfit.fits', oldbg500, oldbg500hdr
fits_read, prefix+'/sp350bg_'+prefix+'_iter'+iterstr+'_gaussfit.fits', oldbg350, oldbg350hdr
fits_read, prefix+'/sp250bg_'+prefix+'_iter'+iterstr+'_gaussfit.fits', oldbg250, oldbg250hdr
fits_read, prefix+'/pa170bg_'+prefix+'_iter'+iterstr+'_gaussfit.fits', oldbg170, oldbg170hdr
fits_read, prefix+'/pa70bg_'+prefix+'_iter'+iterstr+'_gaussfit.fits', oldbg70, oldbg70hdr


sz = size(sp350)
hcongrid, oldbg500, oldbg500hdr, bg500, bg500hdr, /half, outsize=[sz(1),sz(2)]
hcongrid, oldbg350, oldbg350hdr, bg350, bg350hdr, /half, outsize=[sz(1),sz(2)]
hcongrid, oldbg250, oldbg250hdr, bg250, bg250hdr, /half, outsize=[sz(1),sz(2)]
hcongrid, oldbg170, oldbg170hdr, bg170, bg170hdr, /half, outsize=[sz(1),sz(2)]
hcongrid, oldbg70, oldbg70hdr, bg70, bg70hdr, /half, outsize=[sz(1),sz(2)]


;Data - BG, masked to bolocat points (outside that = 0)
fits_read, prefix+'/sp500diff_'+prefix+'_iter'+iterstr+'_gaussfit_conv25.fits',$
  sp500diffmask, sp500diffmaskhdr
fits_read, prefix+'/sp350diff_'+prefix+'_iter'+iterstr+'_gaussfit_conv25.fits',$
  sp350diffmask, sp350diffmaskhdr
fits_read, prefix+'/sp250diff_'+prefix+'_iter'+iterstr+'_gaussfit_conv25.fits',$
  sp250diffmask, sp250diffmaskhdr
fits_read, prefix+'/pa170diff_'+prefix+'_iter'+iterstr+'_gaussfit_conv25.fits',$
  pa170diffmask, pa170diffmaskhdr
fits_read, prefix+'/pa70diff_'+prefix+'_iter'+iterstr+'_gaussfit_conv25.fits',$
  pa70diffmask, pa70diffmaskhdr


npix = n_elements(sp350)
sz = size(sp350)

;get 36" resolution background files and hcongrid to 25" grid
; (for background filler values)
fits_read, prefix+'/bg_temp_'+prefix+'_iter'+iterstr+'_gaussfit_iter_conv36.fits', tempbg, tempbghdr
fits_read, prefix+'/bg_column_'+prefix+'_iter'+iterstr+'_gaussfit_iter_conv36.fits', columnbg, columnbghdr
fits_read, prefix+'/bg_beta_'+prefix+'_iter'+iterstr+'_gaussfit_iter_conv36.fits', betabg, betabghdr
fits_read, prefix+'/bg_temperr_'+prefix+'_iter'+iterstr+'_gaussfit_iter_conv36.fits', bgtemperr, bgtemperrhdr
fits_read, prefix+'/bg_colerr_'+prefix+'_iter'+iterstr+'_gaussfit_iter_conv36.fits', bgcolerr, bgcolerrhdr
fits_read, prefix+'/bg_chisqr_'+prefix+'_iter'+iterstr+'_gaussfit_iter_conv36.fits', bgchisqr, bgchisqrhdr

hcongrid, tempbg, tempbghdr, /half_half, outsize=[sz(1), sz(2)]
hcongrid, columnbg, columnbghdr, /half_half, outsize=[sz(1), sz(2)]
hcongrid, betabg, betabghdr, /half_half, outsize=[sz(1), sz(2)]

hcongrid, bgtemperr, bgtemperrhdr, /half_half, outsize=[sz(1), sz(2)]
hcongrid, bgcolerr, bgcolerrhdr, /half_half, outsize=[sz(1), sz(2)]
hcongrid, bgchisqr, bgchisqrhdr, /half_half, outsize=[sz(1), sz(2)]


;get 36" resolution temp and column files and hcongrid to 25" grid
; (for input guesses)
fits_read, prefix+'/gaussfit_iter_beta175_temp_'+prefix+'_iter'+iterstr+'_conv36.fits', temp36, temp36hdr
fits_read, prefix+'/gaussfit_iter_beta175_column_'+prefix+'_iter'+iterstr+'_conv36.fits', column36, column36hdr

hcongrid, temp36, temp36hdr, /half_half, outsize=[sz(1), sz(2)]
hcongrid, column36, column36hdr, /half_half, outsize=[sz(1), sz(2)]

wz = where(label eq 0, complement=ws)

c = 2.99792458e8 ;m/s

temp25 = fltarr(sz(1), sz(2))
column25 = fltarr(sz(1), sz(2))
label25 = fltarr(sz(1), sz(2))
;beta25 = fltarr(sz(1), sz(2))
chisqr25 = fltarr(sz(1), sz(2))
degfree25 = fltarr(sz(1), sz(2))
temperr25 = fltarr(sz(1), sz(2))
colerr25 = fltarr(sz(1), sz(2))


;Main image (25" resolution)
for jj = 0, sz(2)-1 do begin
    for ii = 0, sz(1)-1 do begin
        if label[ii,jj] eq 0 then begin
            temp25[ii,jj] = tempbg(ii,jj)
            column25[ii,jj] = columnbg(ii,jj)
            label25[ii,jj] = 0
            chisqr25[ii,jj] = bgchisqr(ii,jj)
            temperr25[ii,jj] = bgtemperr(ii,jj)
            colerr25[ii,jj] = bgcolerr(ii,jj)
        endif else begin
            if (sp350diffmask(ii,jj) le 0.) or (sp250diffmask(ii,jj) le 0.) or $
              (pa170diffmask(ii,jj) le 0.) or   finite(pa170diffmask(ii,jj), /nan) or $
              finite(sp250diffmask(ii,jj), /nan) or finite(sp350diffmask(ii,jj), /nan) or $
              finite(sp500diffmask(ii,jj), /nan) then begin 
                temp25[ii,jj] = tempbg(ii,jj)
                column25[ii,jj] = columnbg(ii,jj)
                label25[ii,jj] = 0
                chisqr25[ii,jj] = bgchisqr(ii,jj)
                temperr25[ii,jj] = bgtemperr(ii,jj)
                colerr25[ii,jj] = bgcolerr(ii,jj)
            endif else begin
                inu = [sp500diffmask(ii,jj),sp350diffmask(ii,jj),$
                       sp250diffmask(ii,jj),pa170diffmask(ii,jj),pa70diffmask(ii,jj)]
                err = inu*0.2
 
                inu = inu[nb:nn]
                err = err[nb:nn]
                nun = nu[nb:nn]
                ;set guesses to 36" resolution values
                parinfo_fix2[*].value = [temp36(ii,jj), column36(ii,jj)]
                pfit = mpfitfun('greybody_beta175', nun, inu, err, $
                                /quiet, parinfo=parinfo_fix2, bestnorm=chi, dof=freedeg, covar=comatrix)
                
                chisqr25[ii,jj] = chi
                degfree25[ii,jj] = freedeg
                temperr25[ii,jj] = sqrt(comatrix(0,0))
                colerr25[ii,jj] = sqrt(comatrix(1,1))

                temp25[ii,jj] = pfit(0)
                column25[ii,jj] = pfit(1)
                label25[ii,jj] = 1
            endelse
        endelse
    endfor
endfor

wnan = where( finite(sp500, /nan) or finite(sp350, /nan) or finite(sp250, /nan) or $
              finite(pa170, /nan) )

temp25[wnan] = !values.f_nan
column25[wnan]= !values.f_nan
label25[wnan] = !values.f_nan
temperr25[wnan] = !values.f_nan
colerr25[wnan] = !values.f_nan
chisqr25[wnan] = !values.f_nan


fits_write, prefix+'/gaussfit_iter_beta175_temp_'+prefix+'_iter'+iterstr+'_conv25.fits', temp25, sp500hdr
fits_write, prefix+'/gaussfit_iter_beta175_column_'+prefix+'_iter'+iterstr+'_conv25.fits', column25, sp500hdr
fits_write, prefix+'/gaussfit_iter_beta175_column_'+prefix+'_iter'+iterstr+'_conv25_label.fits', label25, sp500hdr
fits_write, prefix+'/gaussfit_iter_beta175_temperr_'+prefix+'_iter'+iterstr+'_conv25.fits', temperr25, sp500hdr
fits_write, prefix+'/gaussfit_iter_beta175_colerr_'+prefix+'_iter'+iterstr+'_conv25.fits', colerr25, sp500hdr
fits_write, prefix+'/gaussfit_iter_beta175_chisqr_'+prefix+'_iter'+iterstr+'_conv25.fits', chisqr25, sp500hdr


endif



        
END
