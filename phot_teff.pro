;+
; NAME:
;	phot_teff
;
; PURPOSE:
;     >>>>>>>>>>>> turn photometry into Teff and errors<<<<<<<<<<<<<<
;
; CALLING SEQUENCE:
; 
;	phot_Teff,teff,teff_err,vj=vj,vi=vi,rz=rz,rj=rj,jh=jh,e_vj=e_vj,e_vi=e_vi,e_rz=e_rz,e_rj=e_rj,e_jh=e_jh
;
; INPUTS:
;     vj = V-J color (in magnitudes)
;     vi = V-Ic color (I = cousins I)
;     rz = V-J color
;     jh = V-J color
;     e_xx = error on xx color (in magnitudes)
;
; OUTPUTS:
;     teff = effective temperature in kelvin
;     teff_err = error on teff in kelvin
;
;
; RESTRICTIONS:
;     V = johnson V
;     r = SDSS r (apass is ok, but note errors are probably underestimated)
;     z = SDSS z (apass is not ok, color terms)
;     i in vi = cousins I (probably don't use this)
;


;; this is a code to test and make sure it works before I throw it on
;; github. It just pulls down the synthetic colors from the 2015 paper
;; If you feel like running it you need to make sure you have queryvizier
PRO test_photteff

  plotsym,0,/fill
  struct = queryvizier('J/ApJ/804/64/stars',[0,0],1d5,/all) ;; download full catalog
  ;;struct = mrdfits('data.fit',1) ;;
  vj = (struct.VMAG-struct.jmag)
  e_vj = sqrt(struct.e_vmag^2.+struct.e_jmag^2.)
  ;;vi = struct.VMAG-struct.imag
  ;;e_vi = sqrt(struct.e_vmag^2.+struct.e_imag^2.)
  rz = struct.rmag_2-struct.zmag ;; note struct.rmag is Cousins R, not SDSS
  e_rz = sqrt(struct.e_rmag^2.+struct.e_zmag^2.)
  rj = struct.RMAG_2-struct.jmag
  e_rj = sqrt(struct.e_rmag^2.+struct.e_jmag^2.)
  jh = struct.JMAG-struct.HMAG
  e_jh = sqrt(struct.e_jmag^2.+struct.e_hmag^2.)

  big_Teff = dblarr(n_elements(struct))
  big_teff_Err = big_Teff
  for i = 0,n_Elements(struct)-1 do begin
     ;;phot_teff,teff,teff_err,vj=vj[i],rz=rz[i],rj=rj[i],jh=jh[i],e_vj=e_vj[i],e_rz=e_rz[i],e_rj=e_rj[i],e_jh=e_jh[i]
     phot_teff,teff,teff_err,vj=vj[i],rz=rz[i],rj=rj[i],e_vj=e_vj[i],e_rz=e_rz[i],e_rj=e_rj[i]
     big_Teff[i] = teff
     big_teff_err[i] = teff_err
  endfor

  ploterror,struct.teff,big_Teff,struct.e_teff,big_Teff_err,xthick=4,ythick=4,errthick=4,xtitle='Spectroscopic T!Leff!N',ytitle='Photometric T!Leff!N',charsize=1.5,charthick=2,psym=8
  oplot,[0,1d5],[0,1d5],thick=4,color=cgcolor('red'),linestyle=2
  chisq = total((big_teff-struct.teff)^2./(big_Teff_err^2.+struct.e_teff^2.))
  rchisq = chisq/n_elements(big_teff)
  print,'Reduced chi^2: ',string(rchisq,format="(D6.2)")+' (rchisq should be <1 because errors include systematics'
  ;; hopefully this also makes it clear that photometric temperatures
  ;; are quite reliable compared to spectroscopic ones, provided you
  ;; several colors and lots of wavelength baseline (e.g., optical-infrared).
  stop
  
END


PRO phot_Teff,teff,teff_err,vj=vj,vi=vi,rz=rz,rj=rj,jh=jh,e_vj=e_vj,e_vi=e_vi,e_rz=e_rz,e_rj=e_rj,e_jh=e_jh

  nmonte = 1d5
  teff = dblarr(1)
  teff_err = teff
  if n_elements(jh) gt 0 then jh_base = jh

  if n_elements(vj) ne 0 then begin
     col = vj
     if n_elements(e_vj) gt 0 then col_err = e_vj
     if n_elements(jh) gt 0 then jh = jh_base
     if n_elements(col) gt 0 then begin
        if n_elements(jh) gt 0 then begin
           teff = [teff,2.769 -1.421*col+0.4284*col^2.-0.06133*col^3.+0.003310*col^4.+0.1333*jh+0.05416*jh^2]
           if n_elements(e_jh) gt 0 and n_elements(col_err) gt 0 then begin
              col = col+col_err*randomn(seed,nmonte)
              jh = jh+e_jh*randomn(seed,nmonte)
              tmp = 2.769 -1.421*col+0.4284*col^2.-0.06133*col^3.+0.003310*col^4.+0.1333*jh+0.05416*jh^2.
              teff_err = [teff_err,sqrt(robust_sigma(tmp[where(tmp gt 0 and finite(tmp) eq 1)])^2.+(48./3500.)^2.)]
           endif else begin
              teff_err = [teff_err,-99+dblarr(n_elements(jh))]
           endelse
        endif else begin
           teff = [teff,2.840-1.3453*col+0.3906*col^2.-0.0546*col^3.+0.002913*col^4.]
           if n_elements(col_err) gt 0 then begin
              col = col+col_err*randomn(seed,nmonte)
              tmp = 2.840-1.3453*col+0.3906*col^2.-0.0546*col^3.+0.002913*col^4.
              teff_err = [teff_err,sqrt(robust_sigma(tmp[where(tmp gt 0 and finite(tmp) eq 1)])^2.+(55./3500.)^2.)]
           endif else begin
              teff_err = [teff_err,-99+dblarr(n_elements(jh))]
           endelse
        endelse
     endif
  endif

  ;; note this relation uses V-I_cousins, not V-i'
  ;; (sdss). Probably just don't use this because who has cousins I anymore?
  if n_elements(vi) ne 0 then begin
     col = vi
     if n_elements(e_vi) gt 0 then col_err = e_vi
     if n_elements(jh) gt 0 then jh = jh_base
     if n_elements(col) gt 0 then begin
        if n_elements(jh) gt 0 then begin
           teff = [teff,1.568-0.4381*col+0.07749*col^2.-0.005610*col^3.+0.2441*jh-0.09257*jh^2]
           if n_elements(e_jh) gt 0 and n_elements(col_err) gt 0 then begin
              col = col+col_err*randomn(seed,nmonte)
              jh = jh+e_jh*randomn(seed,nmonte)
              tmp = 1.568 -0.4381*col+0.07749*col^2.-0.005610*col^3.+0.2441*jh-0.09257*jh^2
              teff_err = [teff_err,sqrt(robust_sigma(tmp[where(tmp gt 0 and finite(tmp) eq 1)])^2.+(52./3500.)^2.)]
           endif else begin
              teff_err = [teff_err,-99+dblarr(n_elements(jh))]
           endelse
        endif else begin
           teff = [teff,2.455-1.5701*col+0.6891*col^2.-0.1500*col^3.+0.01254*col^4.]
           if n_elements(col_err) gt 0 then begin
              col = col+col_err*randomn(seed,nmonte)
              tmp = 2.455-1.5701*col+0.6891*col^2.-0.1500*col^3.+0.01254*col^4.
              teff_err = [teff_err,sqrt(robust_sigma(tmp[where(tmp gt 0 and finite(tmp) eq 1)])^2.+(53./3500)^2.)]
           endif else begin
              teff_err = [teff_err,-99+dblarr(n_elements(jh))]
           endelse
        endelse
     endif
  endif

  if n_elements(rz) ne 0 then begin
     col = rz
     if n_elements(e_rz) gt 0 then col_err = e_rz
     if n_elements(jh) gt 0 then jh = jh_base
     if n_elements(col) gt 0 then begin
        if n_elements(jh) gt 0 then begin
           teff = [teff,1.384-0.6132*col+0.3110*col^2.-0.08574*col^3.+0.008895*col^4.+0.1865*jh-0.02039*jh^2]
           if n_elements(e_jh) gt 0 and n_elements(col_err) gt 0 then begin
              col = col+col_err*randomn(seed,nmonte)
              jh = jh+e_jh*randomn(seed,nmonte)
              tmp =1.384-0.6132*col+0.3110*col^2.-0.08574*col^3.+0.008895*col^4.+0.1865*jh-0.02039*jh^2
              teff_err = [teff_err,sqrt(robust_sigma(tmp[where(tmp gt 0 and finite(tmp) eq 1)])^2.+(55./3500.)^2.)]
           endif else begin
              teff_err = [teff_err,-99+dblarr(n_elements(jh))]
           endelse
        endif else begin
           teff = [teff,1.547-0.7053*col+0.3656*col^2.-0.1008*col^3.+0.01046*col^4.]
           if n_elements(col_err) gt 0 then begin
              col = col+col_err*randomn(seed,nmonte)
              tmp = 1.547-0.7053*col+0.3656*col^2.-0.1008*col^3.+0.01046*col^4.
              teff_err = [teff_err,sqrt(robust_sigma(tmp[where(tmp gt 0 and finite(tmp) eq 1)])^2.+(58./3500.)^2.)]
           endif else begin
              teff_err = [teff_err,-99+dblarr(n_elements(jh))]
           endelse
        endelse
     endif
  endif

  if n_elements(rj) ne 0 then begin
     col = rj
     if n_elements(e_vi) gt 0 then col_err = e_vi
     if n_elements(jh) gt 0 then jh = jh_base
     if n_elements(col) gt 0 then begin
        if n_elements(jh) gt 0 then begin
           teff = [teff,2.151-1.092*col+0.3767*col^2.-0.06292*col^3.+0.003950*col^4.+0.1697*jh+0.03106*jh^2]
           if n_elements(e_jh) gt 0 and n_elements(col_err) gt 0 then begin
              col = col+col_err*randomn(seed,nmonte)
              jh = jh+e_jh*randomn(seed,nmonte)
              tmp = 2.151-1.092*col+0.3767*col^2.-0.06292*col^3.+0.003950*col^4.+0.1697*jh+0.03106*jh^2
              teff_err = [teff_err,sqrt(robust_sigma(tmp[where(tmp gt 0 and finite(tmp) eq 1)])^2.+(52./3500.)^2.)]
           endif else begin
              teff_err = [teff_err,-99+dblarr(n_elements(jh))]
           endelse
        endif else begin
           teff = [teff,2.445-1.2578*col+0.4340*col^2.-0.0720*col^3.+0.004502*col^4.]
           if n_elements(col_err) gt 0 then begin
              col = col+col_err*randomn(seed,nmonte)
              tmp = 2.445-1.2578*col+0.4340*col^2.-0.0720*col^3.+0.004502*col^4.
              teff_err = [teff_err,sqrt(robust_sigma(tmp[where(tmp gt 0 and finite(tmp) eq 1)])^2.+(58./3500.)^2.)]
           endif else begin
              teff_err = [teff_err,-99+dblarr(n_elements(jh))]
           endelse
        endelse
     endif
  endif
  if n_elements(jh) gt 0 then jh = jh_base
  
  shrink,teff
  shrink,teff_err

  teff*=3500d0
  teff_err*=3500d0

  if n_elements(teff) gt 1 then begin
     teff = wmean(teff,teff_err,err=newteff_Err)
     teff_err = sqrt(newteff_err^2.+60.^2.)
  endif

END
