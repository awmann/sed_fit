
PRO get_flux,struct,system,band,lambda,spec,err,flux,flux_err,fwhm
  COMPILE_OPT idl2, HIDDEN
  if system eq 'Eggen' and band eq 'V' then system = 'Johnson'
  if system eq 'Cousins' and band eq 'U' then system = 'Johnson'
  if system eq 'Cousins' and band eq 'B' then system = 'Johnson'
  if system eq 'Cousins' and band eq 'V' then system = 'Johnson'
  if system eq 'Straizys' then system = 'Vilnius'
  if system eq 'Cape' and (band eq 'V' or band eq 'B') then system = 'Johnson'
  if system eq 'Johnson' and band eq 'I' then begin ;; this system doesnt work
     flux = -99
     flux_err = -99
  endif else begin
     nmonte = 200
     l = where(STRLOWCASE(struct.system) eq STRLOWCASE(system) and STRLOWCASE(struct.band) eq STRLOWCASE(band))
     if l[0] eq -1 then begin
        ;;print,'Warning, missing '+system + ' '+band
        fwhm = 0
        flux = -10
        flux_err = 99
     endif else begin
        flambda = *struct[l].lambda
        ftrans = *struct[l].trans
        gg = where(flambda gt 1d-5 and ftrans gt 1d-5)
        flambda = flambda[gg]
        ftrans = ftrans[gg]
        sp = interpol(spec,lambda,flambda)
        er = interpol(err,lambda,flambda)
        weff = integral(flambda,ftrans)/max(ftrans)
        fwhm = weff/2.
        if system eq 'Gaia' then begin
           vega = mrdfits('alpha_lyr_mod_002.fits',1,header,/silent)
           gvega = 0.023
           vlambda = vega.wavelength ;;/1d4
           vflux = vega.flux
           a = integral(flambda,sp*ftrans*flambda)
           tmp = interpol(vflux,vlambda/1d4,flambda)
           b = integral(flambda,tmp*ftrans*flambda)
           mag = -2.5*alog10(a)+2.5*alog10(b)+gvega
           flux = struct[l].zp*10d0^(-0.4*mag)
           tester = dblarr(nmonte)
           for i = 0,nmonte-1 do begin 
              sp_t = sp+er*randomn(seed,n_elements(sp))
              a = integral(flambda,sp_t*ftrans*flambda)
              mag = -2.5*alog10(a)+2.5*alog10(b)+gvega
              tester[i] = struct[l].zp*10d0^(-0.4*mag)
           endfor
        endif else begin
           flux = integral(flambda,sp*ftrans*flambda)/integral(flambda,ftrans*flambda)
           tester = dblarr(nmonte)
           for i = 0,nmonte-1 do begin 
              sp_t = sp+er*randomn(seed,n_elements(sp)) 
              tester[i] = integral(flambda,sp_t*ftrans)/integral(flambda,ftrans)
           endfor
           ;; is Weff = FWHM?
           ;;result = gaussfit(flambda,ftrans,A,nterms=3)
           ;;fwhm = a[2]
        endelse
        flux_err = stdev(tester)
        if system eq 'Cousins' and band eq 'I' then flux/=0.97 ;; this might be backwards, I think it's right.
     endelse
  endelse

END
