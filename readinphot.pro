;; read in photometry from relevant file
PRO readinphot,r,lambda,spec,err,phot_systems,phot_system,phot_band,phfluxes,phfluxes_err,corr,corr_err,fluxes,fluxes_err,lams,fwhms,mags=mags
  close,/all
  COMPILE_OPT idl2, HIDDEN
  name = r.cns3
  strreplace,name,' ',''
  name = r.cns3
  strreplace,name,' ',''
  strreplace,name,' ',''
  strreplace,name,' ',''
  strreplace,name,' ',''
  strreplace,name,' ',''
  strreplace,name,' ',''
  if strtrim(name,2) eq '' then name = r.name
  file = 'photometry/'+name+'.txt'
  openr,25,file
  line = ''
  lines = file_lines(file)
  phot_system = strarr(1)
  phot_band = strarr(1)
  phfluxes = dblarr(1)
  phfluxes_err = dblarr(1)
  fluxes = dblarr(1)
  fluxes_err = dblarr(1)
  corr = dblarr(1)
  corr_err = dblarr(1)
  lams = dblarr(1)
  fwhms = dblarr(1)
  mags = dblarr(1)
  mags_err = mags
  if lines lt 3 then begin ;; not enough photometry
     print,'there is not enough photometry, the code will certainly crash'
  endif
  for jj = 0,lines-1 do begin
     readf,25,line
     if strmid(line,0,1) ne '#' and strpos(line,' beta') eq -1 and strtrim(line,2) ne '' then begin ;; commented out
        tmp = strsplit(line,' ',/extract)
        phot_err = tmp[n_elements(tmp)-4]
        phot = tmp[n_elements(tmp)-5]
        if finite(phot) eq 1 and (strtrim(tmp[n_elements(tmp)-2],2) ne 'Johnson' or strtrim(tmp[n_elements(tmp)-1],2) ne 'I') then begin
           phot_system = [phot_system,tmp[n_elements(tmp)-2]]  ;; e.g., Johnson
           phot_band = [phot_band,tmp[n_elements(tmp)-1]]      ;; e.g., V

           ;; convert the photometry to a flux here
           ii = n_elements(phot_system)-1
           mags = [mags,phot]
           mags_err = [mags_err,phot_err]
           convert_to_fluxes,phot_systems,phot_system[ii],phot_band[ii],phot,phot_err,phflux,phflux_err,leff
           phfluxes = [phfluxes,phflux]
           phfluxes_err = [phfluxes_err,phflux_err]

           get_flux,phot_systems,phot_system[ii],phot_band[ii],lambda,spec,err,flux,flux_err,fwhm
           if leff gt 2.5 then flux_err = flux*0.02
           fluxes = [fluxes,flux]
           fluxes_err = [fluxes_err,flux_Err]
           fwhms = [fwhms,fwhm]

           tot_err = (phflux/flux)*sqrt((phflux_err/phflux)^2.0+(flux_err/flux)^2.0)
           corr = [corr,phflux/flux]
           corr_err = [corr_err,tot_err]
           lams = [lams,leff]
           ;if strpos(line,'24') ne -1 and strpos(line,'Spitzer') ne -1 then begin
           ;   print,phot,phot_err,phflux,phflux_err,flux,flux_err
           ;   if finite(flux_err) eq 0 then stop
           ;endif
        endif
     endif
  endfor
  ;;
  shrink,corr
  shrink,corr_err
  shrink,phot_band
  shrink,phot_system
  shrink,phot
  shrink,phot_err
  shrink,phfluxes
  shrink,phfluxes_err
  shrink,fluxes
  shrink,fluxes_err
  shrink,lams
  shrink,fwhms
  shrink,mags
  shrink,mags_err

  f = sort(lams)
  corr = corr[f]
  corr_err = corr_err[f]
  phot_band = phot_band[f]
  phot_system = phot_system[f]
  phot = phot[f]
  phot_err = phot_err[f]
  phfluxes = phfluxes[f]
  phfluxes_err = phfluxes_err[f]
  fluxes = fluxes[f]
  fluxes_err = fluxes_err[f]
  lams = lams[f]
  fwhms = fwhms[f]
  close,/all

  ;;forprint,phot_system,phot_band,phfluxes,fluxes
  ;;stop


END
