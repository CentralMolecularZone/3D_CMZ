PRO gauss_iter_wrapper_fixedval, origfile, bgfile, prefix, b, cv, wnan, limit

;prefix = 'l000' 
;limit = 112.

;cv = '36'
;origfile = prefix+'/sp500_'+prefix+'_conv'+cvstr+'.fits'
;bgfile = prefix+'/bgimage0_conv720_'+prefix+'_gaussfit.fits'
;restore, prefix+'/b_'+prefix+'.sav'

cvstr = cv
cvfloat = float(cvstr)

fits_read, origfile, ogdata, oghdr
fits_read, bgfile, bg, bghdr

sz = size(ogdata)
xx = indgen(sz(2))

yvals = fltarr(sz(2))
for jj = 0, sz(2)-1 do yvals[jj] = median(bg(*,jj))

;gauss_fit
gaussfit = fltarr(sz(1), sz(2))
chisqrgaussall = fltarr(sz(1)) ;reduced chisqr

;params3 = [300., 0.25, 50, 0., 300., 50., 0.25]
params2 = [350, 50, -0.1, 0.25]
parinfo = replicate({value:0.D, limited:[0,0], limits:[0.D,0]},7)

for ii = 0, sz(1)-1 do begin
    yy = reform(bg[ii,*])
;    yy = yy(where(finite(yy)))
;    b = b(where(finite(yy)))
    err = yy*0.3
    ygfit = mpfitfun('onedgaussian', b, yy, err, params2, $
                     /quiet, bestnorm=chisqr, dof=dof, /nan)
    chisqrgaussall[ii] = chisqr/dof
    gaussfit[ii,*] = onedgaussian(b, ygfit)
endfor

ygaussfit = mpfitfun('onedgaussian', b, yvals, err, params2, $
                        /quiet, bestnorm=chisqrexpgauss, dof=dofexpgauss)
;gaussfit_med = replicate(1,1014) # onedgaussian(b, ygaussfit) 
;why 1,1014 ?? do sz(1) (x-axis) instead 
gaussfit_med = replicate(1,sz(1)) # onedgaussian(b, ygaussfit) 

fits_write, prefix+'/sp500_'+prefix+'_gaussfit_iterfixed.fits', gaussfit, oghdr
fits_write, prefix+'/sp500_'+prefix+'_gaussfit_median_iterfixed.fits', gaussfit_med, oghdr

gfitdiff = ogdata - gaussfit
gfitdiff_med = ogdata - gaussfit_med

fits_write, prefix+'/gfitdiff_'+prefix+'_iterfixed.fits', gfitdiff, oghdr
fits_write, prefix+'/gfitdiff_med_'+prefix+'_iterfixed.fits', gfitdiff_med, oghdr



;Mask 
label = ogdata
invertlabel = ogdata

; binsize=1
; plothist, gfitdiff, xhist, yhist, bin=binsize, /noplot, /nan
; ;mirror data left of peak
; wp = where(yhist eq max(yhist))
; xwp =xhist(wp)
; xwp = xwp(0)
; datalp = gfitdiff(where(gfitdiff lt (xwp+binsize/2.)))
; plothist, datalp, xhistlp, yhistlp, bin=binsize, /noplot
; nlp = n_elements(xhistlp)
; ;xhistlp = xhistlp[0:nlp-4] ;to get rid of dip in middle ;very small difference
; ;yhistlp = yhistlp[0:nlp-4] ;to get rid of dip in middle
; nlp = n_elements(xhistlp) 
; xhistnew = fltarr(nlp*2)
; yhistnew = fltarr(nlp*2)
; xhistnew[0:nlp-1] = xhistlp
; yhistnew[0:nlp-1] = yhistlp
; revxhist = xhistlp - min(xhistlp) + max(xhistlp)
; revyhist = reverse(yhistlp)
; xhistnew[nlp:(2*nlp -1)] = revxhist
; yhistnew[nlp:(2*nlp -1)] = revyhist

; ;percent of data points that are below the peak
; bgpercent = ( float(n_elements(datalp)) / float(n_elements(gfitdiff)) ) * 100.

; ;fit gaussian to mirrored data
; ;qq = where( (xhistnew gt -100) and (xhistnew lt 100) )
; ;xhistnew = xhistnew[qq]
; ;yhistnew = yhistnew[qq]

; params2 = [max(yhistnew), 0.05, xwp, stddev(gfitdiff, /nan)/4.]
; ;err = yhistnew*0.2
; err = sqrt(yhistnew)
; qq = where(err eq 0)
; if qq(0) ne -1 then err[qq] = 0.1
; allgfit = mpfitfun('onedgaussian', xhistnew, yhistnew, err, params2, $
;                    /quiet)

; sigma = allgfit(3)

; ;limit = xwp + 4.*sigma ;try 4.5 instead of 4 ?? CRAZY
; limit = xwp + float(limitnum)*abs(sigma)
; allgfit(3) = abs(sigma)
; print, limit

limit = limit
print, 'limit: ', limit


;apply limit!
ws = where(gfitdiff ge limit, complement=wn)
label[ws] = 1
label[wn] = 0
invertlabel[ws] = 0
invertlabel[wn] = 1
fits_write, prefix+'/gfit_label_iterfixed_'+prefix+'.fits', label, oghdr

bgimage = double(ogdata*invertlabel)
wz = where(bgimage eq 0)
bgimagenan = bgimage
bgimagenan[wz] = !values.f_nan

fwhm_have = float(cv)
bgimage_conv = convolve_res(720., fwhm_have, bgimagenan, oghdr)
bgimage_conv[wnan] = !values.f_nan
fits_write, prefix+'/bgimagefixed_conv720_'+prefix+'_gaussfit.fits', $
  bgimage_conv, oghdr


END
