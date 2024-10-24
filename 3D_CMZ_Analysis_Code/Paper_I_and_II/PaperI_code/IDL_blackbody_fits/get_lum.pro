

;accidentally deleted previous version (who does that??)
; redone

cv = '36'
fileintemp = 'gc/gaussfit_iter_beta175_temp_gc_itervar_conv'+cv+'.fits'
fileincol = 'gc/gaussfit_iter_beta175_column_gc_itervar_conv'+cv+'.fits'
fileout = 'gc/luminosity_gc_conv'+cv+'.fits'

fits_read, fileintemp, temp, hdr
fits_read, fileincol, col, hdr

;get rid of 100's, make for a weird feature
temp[where(temp eq 100)] = !values.f_nan

;put col in units of m^-2
col = col*1d26

beta = 1.75d

muh2 = 2.8d ;kauffman 2008
mh = 1.673532499d-27 ;kg
h = 6.62606876d-34 ; Js
c = 2.99792458d8 ; m/s
k = 1.3806503d-23 ; J/K

nel = 1d3
min_nu = 1d2
max_nu = 500d10 ;big enough range to be successful, any bigger and it crashes

nu = fill_array(nel, min_nu, max_nu)
dnu = (max_nu-min_nu)/nel

sz = size(temp)
nx = sz(1) & ny = sz(2)

kappanu = dblarr(nel)
lum_per_area = dblarr(nx, ny)
Bnu = dblarr(nel)
Snu = dblarr(nel)
tau_nu = dblarr(nel)


;do this outside the mapping loop
for kk = 0, nel-1 do kappanu[kk] = 0.00200837d * (nu(kk)/340.19d9)^beta ;m^2/kg with gas/dust = 100
;fit using "kappa_fits.pro" with beta = 1.75

for ii = 0, nx-1 do begin
	for jj = 0, ny-1 do begin

for kk = 0, nel-1 do begin
	
	Bnu[kk] = (2d*h*nu(kk)^3d/(c^2d))* (1d/(exp((h*nu(kk))/(k*temp(ii,jj)))-1d)) ;mks
	
	tau_nu[kk] = col(ii,jj) * muh2 * mh * kappanu(kk)
	Snu[kk] = Bnu(kk) * (1 - exp(-1d*tau_nu(kk)))
	lumv = Snu(kk)*dnu 

	lum_per_area[ii,jj] = lum_per_area(ii,jj) + lumv 

endfor
endfor
endfor

pixx = abs(sxpar(hdr, 'CDELT1'))*3600. ;arcseconds per pixel
pixy = abs(sxpar(hdr, 'CDELT2'))*3600. ;arcseconds per pixel

distance = 8.5e3 ;pc, galactic center 

lenx = pixx*distance*3.08e16/206265. ;m
leny = pixy*distance*3.08e16/206265. ;m

area = lenx*leny ;m^2

l_mks = lum_per_area * area ;mks

l = l_mks / 3.386d26 ;Lsun


newhdr = hdr
sxaddpar, newhdr, 'BUNIT', 'Luminosity (Lsun)'

fits_write, fileout, l, newhdr



END


