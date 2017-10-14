PRO struct_radius,str

  nmonte = 100000.
  c = 299792.458D 

  if str.plx gt 0 and str.plx_error gt 0 then begin
     plx = str.plx
     plx_err = str.plx_error
     k = str.k
     k_err = 0.02
     feh = str.irtffeh
     feh = -0.12
     feh_err = 0.1
     plx = abs(plx+plx_err*randomn(seed,nmonte))
     mk = k-5.0*(alog10(1./plx)-1.)
     r2 = (1.9305-0.3466*mk+0.01647*mk^2.)*(1.+feh*0.04458)
     r2 += r2*0.0289*randomn(seed,nmonte)
     r3 = 1.9515-0.3520*mk+0.01680*mk^2.
     print,'Parallax rad: ',median(R2),stdev(r2)
     mass = 5.85825902d-1+3.87151019d-1*mk-1.21729937d-1*mk^2.+1.05529583d-2*mk^3.-2.72615872d-4*mk^4.
     mass += mass*0.018*randomn(seed,nmonte)
     print,'Parallax Mass: ', median(mass),stdev(mass)
     mk+=0.019 ;; http://www.astro.caltech.edu/~jmc/2mass/v3/transformations/ 2MASS -> CIT
     mass = 10.0^(1d-3*(1.8 + 6.12*mk + 13.205*MK^2. - 6.2315*MK^3. + 0.37529*MK^4.))

     den = mass/r2^3.
     tmp = den[sort(den)]
     print,'Density: ',median(den),tmp[n_Elements(tmp)*0.851]-median(den),median(den)-tmp[n_Elements(tmp)*0.159]
  endif

  stop

  ;; first get [Fe/H], use NIR if available, then use optical, then
  ;; JHK relations
  if str.irtffeh ne 0.0 and finite(str.irtffeh) eq 1 then begin
     tmp = *str.irtfspec
     lam = tmp[*,0]
     sp = tmp[*,1]
     err = tmp[*,2]
     rvcorr,str.uh22numspectype,lam,sp,err,3,rv
     lam = lam/(rv/c+1.)

     am_getmetal,lam,sp,err,str.uh22numspectype,'J',feh1,feh1_err,mh1,mh1_err
     am_getmetal,lam,sp,err,str.uh22numspectype,'H',feh2,feh3_err,mh2,mh2_err
     am_getmetal,lam,sp,err,str.uh22numspectype,'K',feh3,feh2_err,mh3,mh3_err
     print,'[Fe/H] = '+string(feh3,format="(D6.2)")+'+/-'+string(feh3_err,format="(D6.2)")
     print,'[M/H] = '+string(mh3,format="(D6.2)")+'+/-'+string(mh3_err,format="(D6.2)")

     feh = feh3
     feh_err = feh3_err;0.09

  endif else begin
     if str.uh22snr gt 50 then begin
        var_r = 1
        var_b = 1
        get_struct_wave_flux,str,wave_r,flux_r,var=var_r,/restwave
        get_struct_wave_flux,str,wave_b,flux_b,var=var_b,/restwave,/b
        wave = [wave_b,wave_r]
        flux = [flux_b,flux_r]
        err = sqrt([var_b,var_r])
        am_getmetal,wave/1d4,flux,err,str.uh22numspectype,'Vis',feh,feh_err,mh,mh_err
        feh_err = sqrt(feh_err^2.+0.11^2.)
        mh_err = sqrt(mh_err^2.+0.11^2.)
        print,'[Fe/H] = '+string(feh,format="(D6.2)")+'+/-'+string(feh_err,format="(D6.2)")
        print,'[M/H] = '+string(mh,format="(D6.2)")+'+/-'+string(mh_err,format="(D6.2)")
     endif
  endelse
  ;feh = -0.03
  ;feh_err = 0.09
  ;str.uh22model_values[6] = 3307
  
  ;; next get the Teff, ideally we should use the teff from Spectra, 
  ;; but calculate everything you can
  vj = str.v-str.j
  if vj gt 2.5 and vj lt 8 and str.v gt 0 and str.j gt 0 and finite(str.v) eq 1 and finite(str.j) eq 1 then begin
     teff = (2.840-1.3453*vj+0.3906*vj^2.-0.0546*vj^3.+0.002913*vj^4.)*3500.
     teff_err = sqrt(55^2.+60^2.);; not quite right, should mc this
     print,'T_vj = '+string(teff,format="(I4)")+'+/-'+string(teff_err,format="(I4)")
  endif
  rj = str.v-str.j
  if rj gt 2.5 and rj lt 8 and str.r gt 0 and str.j gt 0 and finite(str.r) eq 1 and finite(str.j) eq 1 then begin
     teff = (2.445-1.2578*rj+0.4340*rj^2.-0.0720*rj^3.+0.004502*rj^4.)*3500.
     teff_err = sqrt(58^2.+60^2.);; not quite right, should mc this
     print,'T_rj = '+string(teff,format="(I4)")+'+/-'+string(teff_err,format="(I4)")
  endif
  if str.uh22model_Values[6] gt 20 then begin
     teff = str.uh22model_values[6]
     teff_err = 60.
     print,'T_sp = '+string(teff,format="(I4)")+'+/-'+string(teff_err,format="(I4)")
  endif 
     

  t=(teff/3500.)+(teff_err/3500.)*randomn(seed,nmonte)
  f=feh+feh_err*randomn(seed,nmonte)
  
  r1 = (16.7700-54.3210*t+57.6627*t^2.-19.6994*t^3.)*(1+0.4565*f)
  f_tmp=feh+0.08*randomn(seed,nmonte)
  t_tmp=(teff/3500.)+(60./3500.)*randomn(seed,nmonte)
  r1_old = (16.7700-54.3210*t_tmp+57.6627*t_tmp^2.-19.6994*t_tmp^3.)*(1+0.4565*f_tmp)
  extra = sqrt(robust_sigma(r1[where(finite(r1) eq 1)])^2.-robust_sigma(r1_old[where(finite(r1_old) eq 1)])^2.)
  if extra lt 0 or finite(extra) eq 0 then extra = 0
  r1_err = sqrt((0.093*median(r1))^2.+extra^2.)

  mass = 0.0078368671+0.48699518*r1+1.6777663*r1^2.-1.2430190*r1^3.
  mass_err = sqrt(stdev(mass[where(finite(mass) eq 1)])^2.+(median(mass[where(finite(mass) eq 1)])*0.1)^2.)
  rho = mass/r1^3.
  rho_err = 0.06*3*median(rho[where(finite(rho) eq 1)])
  rho_err = stdev(rho[where(finite(rho) eq 1)])
  
  r2 = (10.5440-33.7546*t+35.1909*t^2-11.5928*t^3.)
  t_tmp=(teff/3500.)+(60./3500.)*randomn(seed,nmonte)
  r2_old = (10.5440-33.7546*t_tmp+35.1909*t_tmp^2-11.5928*t_tmp^3.)
  extra = sqrt(robust_sigma(r2)^2.-robust_sigma(r2_old)^2.)
  if extra lt 0 then extra = 0
  r2_err = sqrt((0.134*median(r2))^2.+extra^2.)

  print,'R* = '+string(median(r1[where(finite(r1) eq 1)]),format="(D6.3)")+'+/-'+string(r1_err,format="(D6.3)")
  print,'M* = '+string(median(mass[where(finite(mass) eq 1)]),format="(D6.3)")+'+/-'+string(mass_err,format="(D6.3)")
  print,'\rho* = '+string(median(rho[where(finite(rho) eq 1)]),format="(D6.3)")+'+/-'+string(rho_err,format="(D6.3)")
  ;;print,'R* =
  ;;'+string(median(r2),format="(D6.3)")+'+/-'+string(r2_err,format="(D6.3)")

  logg = alog10(6.6743d-8*mass*1.989d33/(r1*6.955d10)^2.)
  print,median(logg)
  

END


;; returns a modified wavelength array based on the rv correction
PRO rvcorr,spectype,lambda,spec,error,order,rv

  restore,'~/Dropbox/IRTF_obs/tellcor.dat'
  tellerror = tellcorr.mastererror
  
  c = 299792.458D
  oldlambda = lambda
  oldspec = spec
  olderror = error

  if spectype lt 4.0 then comp = mrdfits('~/Dropbox/Structures/All_fits_091201/M1.5V_HD36395.fits',0,header,/silent)
  if spectype ge 4.0 and spectype lt 7.0 then comp = mrdfits('~/Dropbox/Structures/All_fits_091201/M5V_Gl51.fits',0,header,/silent)
  if spectype ge 7.0 then comp = mrdfits('~/Dropbox/Structures/All_fits_091201/M9V_LHS2065.fits',0,header,/silent)
  oldcomplambda = comp[*,0] 
  oldcompspec = comp[*,1]
  oldcomperror = comp[*,2]

  complambda = oldcomplambda
  compspec = oldcompspec
  comperror = oldcomperror

  ;;plot,lambda,spec
  case order of
     3: good = where(oldlambda gt 2.1 and oldlambda lt 2.3)   ;
     4: good = where(oldlambda gt 1.5 and oldlambda lt 1.73)  ;
     5: good = where(oldlambda gt 1.17 and oldlambda lt 1.31) ; where(tellerror/median(tellerror) lt 1.2)
     6: good = where(oldlambda gt 0.96 and oldlambda lt 1.10) ; good = where(tellerror/median(tellerror) lt 1.2)
     0: good = where(tellerror/median(tellerror) lt 1.2)
  endcase
  
  lambda = lambda[good]
  spec = oldspec[good]
  error = olderror[good]
  compspec = interpol(oldcompspec,oldcomplambda,lambda)

  ;;result = linfit(lambda,spec)
  ;;spec = spec/(result[1]*lambda)
  ;;result = linfit(lambda,compspec)
  ;;compspec = compspec/(result[1]*lambda)

  limits = [median(lambda)-stdev(lambda),median(lambda)+stdev(lambda)]

  compute_xcor_rv,lambda, spec,'M1.5', rv, wavelim=[min(lambda),max(lambda)],norm_reg=limits,twave=lambda, tflux=compspec, maxshift=75, showplot=0
  

  lambda = oldlambda
  spec = oldspec
  error = olderror

END
