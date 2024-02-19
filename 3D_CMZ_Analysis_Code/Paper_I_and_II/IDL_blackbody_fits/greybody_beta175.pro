FUNCTION greybody_beta175, nu, p

;input p = [T, Nh2], units [K, x10^22 cm^-2]

Tdust = p(0)  ;units K
Nh2 = p(1)*1e26    ;units #/m^2
;beta = p(2)   ;unitless
beta = 1.75

;kappanu = 0.04 * (nu/5.936e11)^beta
kappanu = 0.0200837 * (nu/340.19e9)^beta ;with gas/dust = 100
;fit using "kappa_fits.pro" with beta = 1.75
kappanu *= 1/10. ;cm^2/g to m^2/kg

muh2 = 2.8 ;kauffman 2008
mh = 1.673532499d-27 ;kg
h = 6.62606876d-34 ; Js
c = 2.99792458d8 ; m/s
k = 1.3806503d-23 ; J/K

Bnu = (2.*h*nu^3./(c^2d))* (1d/(exp((h*nu)/(k*Tdust))-1d)) ;mks
Bnu *= 1d26 ;convert to Jy

;optically thin dust
;Snu = Bnu * Nh2 * muh2 * mh * kappanu; Jy/Sr
tau_nu = Nh2 * muh2 * mh * kappanu
Snu = Bnu * (1 - exp(-1.*tau_nu))

Snu *= 1d-6 ;to MJy/Sr

;stop
return, Snu

END

