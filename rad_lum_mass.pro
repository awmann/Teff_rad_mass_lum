;+
; NAME:  
;       rad_lum_mass
; PURPOSE:
;        Convert effective temperature and error in effective
;        temperature to radius, mass, and luminosity (with errors)
;        using the empirical formulas from Mann et al. (2013c).
;
;        Program for calculating T_eff from spectra is coming!
;
;        Contact me at amann@astro.as.utexas.edu if you have problems.;        
;
; CALLING SEQUENCE: 
;       rad_mass_lum,teff,teff_err,silent=silent
;
; INPUT PARAMETERS: 
;       teff     Effective temperature in K
;       teff_err Error on effective temperature
;
; OUTPUT PARAMETERS: 
;       The function returns a 2x3 array containing:
;       [[Radius, Radius_err],[Luminosity,Luminosity_err],[Mass,Mass_err]]
;       Output quantities are all in solar units
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;       /silent   supresses some output
;   
; EXAMPLES:
;
;
; DEPENDENCIES:
;
;       
; HISTORY: 
;       Routine written. A. Mann 12/15/2014
;-

FUNCTION rad_lum_mass,teff,teff_err,silent=silent

  if n_elements(silent) eq 0 then silent = 0 else silent = 1

  if n_elements(teff) eq 0 or n_elements(teff_err) eq 0 then begin
     print,'Syntax: rad_lum_mass,teff,teff_err'
     return,-1
  endif

  if teff lt 3180 or teff gt 4773 then begin
     print,'Warning, temperature is outside the range of calibrators!'
     print,'Results are probably not reliable.'
  endif

  readcol,'~/Dropbox/Radii/Relations.txt',/silent,type,a,b,c,d,rchisq,sig,fit_err2,format="a,d,d,d,d,d,d,d"
  readcol,'~/Dropbox/Radii/T_R.dat',newx,t_r_err,/silent
  readcol,'~/Dropbox/Radii/T_L.dat',newx,t_l_err,/silent
  readcol,'~/Dropbox/Radii/T_M.dat',newx,t_m_err,/silent

  values = dblarr(2,3) 
  newx = generatearray(3000,5000,100)
  ;; y is a generic quantity (could be mass, radius, or luminosity)
  for i = 0,n_elements(a)-2 do begin ;; -2 because the last relation is not needed for this program
     y = a[i]+b[i]*teff+c[i]*teff^2.0+d[i]*teff^3.0
     ;; error analysis: the error in y is sqrt((dy/dT)^2.0*\sigmaT^2.0
     ;; + \sigma_Relation^2.0) 
     newy = a[i]+b[i]*newx+c[i]*newx^2.0+d[i]*newx^3.0

     dy = deriv(newx,newy)
     dy_T = interpol(dy,newx,teff)

     case i of 
        0: begin ;; Teff-Radius
           fit_err = interpol(t_r_err,newx,teff)
        end
        1: begin ;; Teff-Luminosity
           fit_err = interpol(t_l_err,newx,teff)
        end
        2: begin ;; Teff-Mass
           fit_err = interpol(t_m_err,newx,teff)
        end
     endcase
     sigma = sqrt(dy_t^2.0*teff_err^2.0 + fit_err^2.0)
     values[*,i] = [y,sigma]
  endfor

  ;; fundamental lower mass error limit of 10%
  if values[1,2]/values[0,2] lt 0.1 then values[1,2] = values[0,2]*0.1 ;; this should rarely execute
  if silent eq 0 then begin
     print,'Radius: '+string(values[0,0],format="(D6.3)")+'+/-'+string(values[1,0],format="(D6.3)")
     print,'Luminosity: '+string(values[0,1],format="(D6.3)")+'+/-'+string(values[1,1],format="(D6.3)")
     print,'Mass: '+string(values[0,2],format="(D6.3)")+'+/-'+string(values[1,2],format="(D6.3)")
  endif
  return,values

END
