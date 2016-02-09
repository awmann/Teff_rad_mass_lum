;+
; NAME:
;	phot_teff
;
; PURPOSE:
;     >>>>>>>>>>>> turn Teff into radius and errors<<<<<<<<<<<<<<
;
; CALLING SEQUENCE:
; 
;	rad_teff,teff,teff_err,rad,rad_err,feh=feh,feh_Err=feh_err
;
; INPUTS:
;     teff = effective temperature
;     teff_err = error on teff
;     feh = metallicity (optional but improves precision)
;     feh_err = error on feh
;
; OUTPUTS:
;     rad = radius in solar units
;     rad_err = error on radius
;
;
; RESTRICTIONS:
;

PRO test_radteff

  struct = queryvizier('J/ApJ/804/64/stars',[0,0],1d5,/all) ;; download full catalog
  plotsym,0,/fill
  big_rad = dblarr(n_elements(struct))
  big_rad_Err = big_rad
  for i = 0,n_Elements(struct)-1 do begin
     ;;rad_Teff,struct[i].teff,struct[i].e_teff,rad,rad_err,feh=struct[i]._FE_H_,e_feh=struct[i].E__FE_H_
     rad_Teff,struct[i].teff,struct[i].e_teff,rad,rad_err
     big_rad[i] = rad
     big_rad_err[i] = rad_err
  endfor

  ploterror,struct.r,big_rad,struct.e_r,big_rad_err,psym=8,/xstyle,/ystyle
  oplot,[0,1d5],[0,1d5],thick=2,color=cgcolor('red'),linestyle=2
  chisq = total((struct.r-big_rad)^2./(big_rad_err^2.+struct.e_r^2.))
  rchisq = chisq/n_elements(struct)
  print,'Reduced Chi^2: '+string(rchisq,format="(D6.2)")
  ;; again this should be <1, because we are actually quite
  ;; conservative with errors
  stop

END

PRO rad_teff,teff,teff_err,rad,rad_err,feh=feh,e_feh=e_feh

  nmonte = 1d5
  if n_elements(teff_err) eq 0 then teff_err = 0d0
  if n_elements(feh) eq 0 then begin
     t = teff/3500d0
     rad = 10.5440-33.7546*t+35.1909*t^2.-11.5928*t^3.
     tmp = teff+randomn(seed,nmonte)*teff_err	
     t = tmp/3500d0				
     rad2 = 10.5440-33.7546*t+35.1909*t^2.-11.5928*t^3.
     rad_err = sqrt(robust_sigma(rad2)^2.+(0.134*rad)^2.)
  endif else begin
     t = teff/3500d0
     rad = (16.7700-54.3210*t+57.6627*t^2.-19.6994*t^3.)*(1d0+0.4565*feh)
     tmp = teff+randomn(seed,nmonte)*teff_err	
     t = tmp/3500d0
     tmp = feh+randomn(seed,nmonte)*e_feh
     rad2 = (16.7700-54.3210*t+57.6627*t^2.-19.6994*t^3.)*(1d0+0.4565*tmp)
     rad_err = sqrt(robust_sigma(rad2)^2.+(0.093*rad)^2.)
  endelse

END
