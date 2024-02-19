

;accidentally deleted previous version (who does that??)
; redone

cv = '25'
fileintemp = 'gc/gaussfit_iter_beta175_temp_gc_itervar_conv'+cv+'.fits'
fileout = 'gc/luminosity_gc_conv'+cv+'.fits'

fits_read, fileintemp, temp, hdr

;get rid of 100's, make for a weird feature
temp[where(temp eq 100)] = !values.f_nan


beta = 1.75

muh2 = 2.8d ;kauffman 2008
mh = 1.673532499d-27 ;kg
h = 6.62606876d-34 ; Js
c = 2.99792458d8 ; m/s
k = 1.3806503d-23 ; J/K


pixx = abs(sxpar(hdr, 'CDELT1'))*3600. ;arcseconds per pixel
pixy = abs(sxpar(hdr, 'CDELT2'))*3600. ;arcseconds per pixel

distance = 8.5e3 ;pc, galactic center 

lenx = pixx*distance*3.08e16/206265. ;m
leny = pixy*distance*3.08e16/206265. ;m

area = lenx*leny ;m^2

;analytic solution to the integral of Bnu over nu
; times area = luminosity (energy s^-2 Sr^-1)
; Page 101 of NB, March 2014
Lprefix = (2.*(!dpi)^4. * k^4.) / (15d*h^3.*c^2.)
L_mks = Lprefix * area * temp^4. ;mks

L = L_mks/3.839d26 ;Lsun


newhdr = hdr
sxaddpar, newhdr, 'BUNIT', 'Luminosity (Lsun)'


fits_write, fileout, L, newhdr




END


