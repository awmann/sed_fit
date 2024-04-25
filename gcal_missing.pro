PRO iras
  rsun = 6.957d10
  lsun = 3.828d33
  ; restore,'young.dat'
  ; restore,'mdallax.dat'
  ; yo = young[where(young.uh22numspectype gt -0.4 and young.uh22numspectype lt 1.5 and young.uh22snr gt 50 and young.irtfsnr gt 50)]
  ; y = yo[0]
  ; m = mdallax[where(abs(mdallax.uh22model_values[6]-4000) lt 200 and mdallax.uh22spec_b ne ptr_new())]
  ; m = m[sort(m.uh22model_values[6])]
  ; m_save = m
  ; young_good = young[where(abs(young.uh22model_values[6]-4000) lt 200 and young.uh22snr gt 50 and young.irtfsnr gt 100 and young.uh22spec_b ne ptr_new())]
  ; young_good = young_good[where(young_good.name ne 'PM_I03185+4432')]
  ; young_good = young_good[sort(young_good.uh22model_values[6])]

  ; save,m,filename='m.idlsav'
  ; save,young_good,filename='young_good.idlsav'
  ; save,y,filename='y.idlsav'
  ; models = mrdfits('Models_CIFIST_Aug2013_filler.fits',1,/silent)
  ; l = where(models.teff gt 3700 and models.teff lt 4500 )
  ; models.spectrum = (models.spectrum)[l,*]
  ; models.TEFF = (models.TEFF)[l]
  ; models.LOGG = (models.LOGG)[l]
  ; models.METAL = (models.METAL)[l]
  ; models.A_FE = (models.A_FE)[l]
  ; A = {header:models.header, Spectrum:(models.spectrum)[l,*], TEFF:(models.TEFF)[l],logg:(models.LOGG)[l],metal:(models.METAL)[l],A_FE:(models.A_FE)[l]}
  ; models = A
  ; save,models,filename='models.idlsav'
  restore,'m.idlsav'
  restore,'young_good.idlsav'
  restore,'y.idlsav'
  restore,'phot_systems.dat'
  print,n_elements(m),n_elements(young_good)

  fbolgood = []
  lgood = []
  rad_Stef_good = []
  teff_good = []
  r_good = []
  ebvgood = []
  teff_fit_good = []
  add = 1 ;; add models (set to 0 to make it run really fast, 1 for accuracy, 0 is good for testing).
  print,'SpT   rchisq  E(B-V)  Fbol  e_Fbol ModT Lum  e_Lum  R_stef e_R_stef T_IR  R_IR  e_R_IR  savemetric radoff'
  plx =  6.2474 	
  plx_err = 0.027
  cut1 = 3
  cut2 = 3
  cuttot = 5
  add = 1
  bestchisq = 3
  for i = 0,n_elements(m)-1 do begin
    ;;print,m[i].name
    for ebv = 0.7,0.96,0.01 do begin
      spt = m[i].uh22spectype
      strreplace,spt,' ',''
      ; y = yo[i]
      struct_assign,m[i],y
      y.name = 'IRAS04125+2902'
      y.othername = 'IRAS04125+2902'
      y.cns3 = 'IRAS04125+2902'
      name = 'IRAS04125+2902'
      y.ra = 063.9282796869600
      y.dec = +29.1666201244000

      y.irtfsnr = 200

      ;; redden the template
      get_Struct_wave_flux,y,wave,flux,var=var
      ccm_unred,wave,flux,-1d0*ebv,redspec
      y.uh22spec = ptr_new(redspec)
      varb = 1
      get_Struct_wave_flux,y,waveb,fluxb,var=varb,/blue
      ccm_unred,wave,fluxb,-1d0*ebv,redspecb
      y.uh22spec_b = ptr_new(redspecb)

      tmp = *y.irtfspec
      sp = tmp[*,1]
      wave = tmp[*,0]
      err = tmp[*,2]
      l = where(wave lt max(wave)-0.1)
      sp = sp[l]
      wave = wave[l]
      err = err[l]
      ccm_unred,wave*1d4,sp,-1d0*ebv,redspecir
      tmp = [[wave],[redspecir],[err]]
      y.irtfspec = ptr_new(tmp)

      gcal,y,phot_systems,finallambda,finalspec,finalerr,ffbol,finalfbol_err,rchisq,rchisq_nowise=rchisq_nowise,redo=redo,add=add,hot=0,PXRANGE=[0.35,10],/quiet,/nostop,/wexcess,ebv=ebv,fitteff=fitteff
      ll = where(finalerr eq median(finalerr))
      ccm_unred,finallambda*1d4,finalspec,ebv,unredspec
      ; unredspec[ll] = finalspec[ll]
      ffbol = integral(finallambda,unredspec)*1d12
      finalfbol_err = finalfbol_err*3.0*(integral(finallambda,unredspec)/integral(finallambda,finalspec))

      nmonte = 1d5
      fbol = ffbol*1d-12
      fbol_err = finalfbol_err*1d-12 + ffbol*1d-12*0.01
      c1 = (4.*!pi * 1d4) / 3.828d33 ;;3.939d33
      c2 = 3.08567758d18
      const = c1*c2^2.0
      dist = 1000d0/(plx+plx_err*randomn(seed,nmonte))
      fbol = fbol+fbol_err*randomn(seed,nmonte)
      luminosity = const*(fbol*(dist)^2.0)
      logluminosity = alog10(luminosity)

      ;; ok, now do the r^2/d^2 trick:

      loggrange = [3.4,4.1]
      metalrange = [-0.1,0.1]
      teffrange = [3750,4100]
      models = mrdfits('~/Dropbox/Radii/Models_CIFIST_Aug2013_filler.fits',1,/silent)

      modelspectra = models.spectrum
      modelheader = models.header
      lambda_m = (sxpar(modelheader,'LAMBDA_0') + findgen(sxpar(modelheader,'NLAMBDA'))*sxpar(modelheader,'D_LAMBDA'))/1d4
      teff = models.teff
      logg = abs(models.logg)
      metal = models.metal
      savemetric = 1d0
      ;;print,'Teff, r, r_err, scale, scatter/scale'
      set_plot,'x'
      !p.multi=[0,1,2]
      matchrange = [1.0,2.4]
      for j = 0,n_elements(models.teff)-1 do begin
        if (models.teff)[j] ge teffrange[0] and (models.teff)[j] le teffrange[1] and $
        (models.logg)[j] ge loggrange[0] and (models.logg)[j] le loggrange[1] and $
        (models.metal)[j] ge metalrange[0] and (models.metal)[j] le metalrange[1] then begin
        sp = modelspectra[j,*]

        l = where(lambda_m gt matchrange[0] and lambda_m lt matchrange[1])
        ll = where(finallambda gt matchrange[0] and finallambda lt matchrange[1])
        lam = finallambda[ll]
        fspec = unredspec[ll]
        tmp = interpol(sp[l],lambda_m[l],lam)
        fac = smooth((fspec/tmp),20)
        scale = median(fac)
        scatter = robust_sigma(fac)
        if scatter le 0 then scatter = stdev(fac)

        scale_tmp = scale + scatter*randomn(seed,nmonte)
        plx_tmp = plx+plx_err*randomn(seed,nmonte)
        d = (1000d0/plx_tmp)*30856776000000d0
        r = sqrt(scale_tmp*d^2d0)/696000d0
        if scatter/scale lt savemetric then begin
          saveteff = (models.teff)[j]
          saver = median(r)
          saver_err = stdev(r)
          savescale = scale
          savemetric = scatter/scale
        endif
        if j gt 0 then begin ;; do half teff
          sp1 = modelspectra[j,*]
          loc = (wherE(models.logg eq (models.logg)[j] and models.metal eq (models.metal)[j] and models.teff eq (models.teff)[j]+100))[0]
          sp2 = modelspectra[loc,*]
          sp = (sp1+sp2)/2d0
          l = where(lambda_m gt matchrange[0] and lambda_m lt matchrange[1])
          ll = where(finallambda gt matchrange[0] and finallambda lt matchrange[1])
          lam = finallambda[ll]
          fspec = unredspec[ll]
          tmp = interpol(sp[l],lambda_m[l],lam)
          fac = smooth((fspec/tmp),20)
          scale = median(fac)
          scatter = robust_sigma(fac)
          if scatter le 0 then scatter = stdev(fac)

          scale_tmp = scale + scatter*randomn(seed,nmonte)
          plx_tmp = plx+plx_err*randomn(seed,nmonte)
          d = (1000d0/plx_tmp)*30856776000000d0
          r = sqrt(scale_tmp*d^2d0)/696000d0
          if scatter/scale lt savemetric then begin
            saveteff = ((models.teff)[j]+(models.teff)[loc])/2.0
            saver = median(r)
            saver_err = stdev(r)
            savescale = scale
            savemetric = scatter/scale
          endif

          ;; and quarter teff
          sp1 = modelspectra[j,*]
          loc = (wherE(models.logg eq (models.logg)[j] and models.metal eq (models.metal)[j] and models.teff eq (models.teff)[j]+100))[0]
          sp2 = modelspectra[loc,*]
          sp = (0.25*sp1+0.75*sp2)
          l = where(lambda_m gt matchrange[0] and lambda_m lt matchrange[1])
          ll = where(finallambda gt matchrange[0] and finallambda lt matchrange[1])
          lam = finallambda[ll]
          fspec = unredspec[ll]
          tmp = interpol(sp[l],lambda_m[l],lam)
          fac = smooth((fspec/tmp),20)
          scale = median(fac)
          scatter = robust_sigma(fac)
          if scatter le 0 then scatter = stdev(fac)

          scale_tmp = scale + scatter*randomn(seed,nmonte)
          plx_tmp = plx+plx_err*randomn(seed,nmonte)
          d = (1000d0/plx_tmp)*30856776000000d0
          r = sqrt(scale_tmp*d^2d0)/696000d0
          if scatter/scale lt savemetric then begin
            saveteff = (0.25*(models.teff)[j]+0.75*(models.teff)[loc])
            saver = median(r)
            saver_err = stdev(r)
            savescale = scale
            savemetric = scatter/scale
          endif

          sp1 = modelspectra[j,*]
          loc = (wherE(models.logg eq (models.logg)[j] and models.metal eq (models.metal)[j] and models.teff eq (models.teff)[j]+100))[0]
          sp2 = modelspectra[loc,*]
          sp = (0.75*sp1+0.25*sp2)
          l = where(lambda_m gt matchrange[0] and lambda_m lt matchrange[1])
          ll = where(finallambda gt matchrange[0] and finallambda lt matchrange[1])
          lam = finallambda[ll]
          fspec = unredspec[ll]
          tmp = interpol(sp[l],lambda_m[l],lam)
          fac = smooth((fspec/tmp),20)
          scale = median(fac)
          scatter = robust_sigma(fac)
          if scatter le 0 then scatter = stdev(fac)

          scale_tmp = scale + scatter*randomn(seed,nmonte)
          plx_tmp = plx+plx_err*randomn(seed,nmonte)
          d = (1000d0/plx_tmp)*30856776000000d0
          r = sqrt(scale_tmp*d^2d0)/696000d0
          if scatter/scale lt savemetric then begin
            saveteff = (0.75*(models.teff)[j]+0.25*(models.teff)[loc])
            saver = median(r)
            saver_err = stdev(r)
            savescale = scale
            savemetric = scatter/scale
          endif

        endif
      endif
    endfor
    
    teff = m[i].uh22model_values[6]+75d0*randomn(seed,nmonte)
    sigma = 5.6704d-5
    rad_stefan = sqrt((luminosity*3.828d33)/(teff^4.*4.*!pi*sigma)) / (6.955d10)

    ;; save only if it's half-decent
    radoff = abs(median(rad_stefan)-saver)/sqrt(saver_err^2.+stdev(rad_stefan)^2.)
    if rchisq_nowise lt cut1 and radoff lt cut2 and (radoff+rchisq_nowise) lt cuttot then begin
      str = 'mv Phot_fits/'+name+'_'+string(add,format="(I1)")+'.eps Phot_fits/'+name+'_'+spt+'_'+strtrim(string(ebv,format="(D6.2)"),2)+'_'+strtrim(string(rchisq_nowise*10,format="(I2)"),2)+'.eps'
      spawn,str
      str2 = 'pstopdf Phot_fits/'+name+'_'+spt+'_'+strtrim(string(ebv,format="(D6.2)"),2)+'_'+strtrim(string(rchisq_nowise*10,format="(I2)"),2)+'.eps Phot_fits/'+name+'_'+spt+'_'+strtrim(string(ebv,format="(D6.2)"),2)+'_'+strtrim(string(rchisq_nowise*10,format="(I2)"),2)+'.pdf'
      spawn,str2
      print,spt,     rchisq,  ebv,    ffbol,  finalfbol_err,m[i].uh22model_values[6],median(luminosity),stdev(luminosity),median(rad_stefan),stdev(rad_stefan),saveteff,saver,saver_err,savemetric,radoff,format=$
      "(A-5,1x,      D5.2,2x, D5.2,2x,D6.3,1x,D6.3,2x,      I4,                      D6.3,1x,D6.3,2x,                     D5.2,1x,D5.2,2x,     I4,2x,D5.2,1x,D5.2,2x,D5.2,3x,D5.3)"
      fbolgood = [fbolgood,ffbol]
      lgood = [lgood,median(luminosity)]
      rad_Stef_good = [rad_Stef_good,median(rad_stefan)]
      teff_good = [teff_good,saveteff]
      r_good = [r_good,saver]
      ebvgood = [ebvgood,ebv]
      teff_fit_good = [teff_fit_good,m[i].uh22model_values[6]]
      if rchisq_nowise lt bestchisq then begin
        bestchisq = rchisq_nowise
        ;; write unreddened spectrum to a file
        tmp = strtrim(string(finallambda*1d4,format="(D9.2)"),2)
        forprint,tmp,unredspec,textout='IRAS_outspec.ascii',comment='Wave Spec'

      endif
    endif

  endfor
endfor

print,'Fbol: '+string(mean(fbolgood),format="(D5.3)")+'+/-'+string(stdev(fbolgood),format='(D5.3)')
print,'Lum: '+string(mean(lgood),format='(D5.3)')+'+/-'+string(stdev(lgood),format='(D5.3)')
print,'T fit: '+string(mean(teff_fit_good),format='(I4)')+'+/-'+string(stdev(teff_fit_good),format='(I4)')
print,'Rad S-B: '+string(mean(rad_Stef_good),format='(D5.3)')+'+/-'+string(stdev(rad_Stef_good),format='(D5.3)')
print,'T IRFM: '+string(mean(teff_good),format='(I4)')+'+/-'+string(stdev(teff_good),format='(I4)')
print,'Rad IRFM: '+string(mean(r_good),format='(D5.3)')+'+/-'+string(stdev(r_good),format='(D5.3)')
print,'AV: '+string(mean(ebvgood*3.1),format='(D5.3)')+'+/-'+string(stdev(ebvgood*3.1),format='(D5.3)')
stop


END



PRO gcal,r,phot_systems,finallambda,finalspec,finalerr,finalfbol,finalfbol_err,rchisq,redo=redo,$
  rchisq_nowise=rchisq_nowise,goodpoints=goodpoints,add=add,quiet=quiet,dup=dup,hot=hot,ebv=ebv,$
  pxrange=pxrange,noplot=noplot,fitteff=fitteff,wexcess=wexcess,nostop=nostop,$
  photometrycolor=photometrycolor,nolegend=nolegend

  set_plot,'x'
  if n_elements(ebv) eq 0 then ebv = 0
  if n_elements(dup) eq 0 then dup = 0
  if n_elements(wexcess) eq 0 then wexcess = 0
  if n_elements(quiet) eq 0 then quiet = 0
  if n_elements(redo) eq 0 then redo = 0
  if n_elements(add) eq 0 then add = 1
  if n_elements(hot) eq 0 then hot = 0
  if n_elements(noplot) eq 0 then noplot = 0
  if n_elements(nostop) eq 0 then nostop = 0
  device,retain=2
  close,/all
  charthick = 3
  charsize = 1.25
  xthick = 5
  ythick = 5
  thick = 3.0
  errthick = 3.0
  plotsym,0,/fill
  xrange = [0.3,4.0]            ;[0.25,0.6];
  xrange_blue = [0.05,0.6]
  xrange_red = [2.0,30]
  !except = 0
  !p.font = 0
  set_plot,'x'
  symsize=1.1

  !p.multi=[0,1,1]
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
  rem = ''

  if quiet ne 1 then print,r.cns3,' ',r.othername,' ',r.name

  create_photometry_file,r,redo,rem    ;; create the fancy photometry file

  satisfied = 0
  counter = 0
  r_fix = 1.0
  b_fix = 1.0
  r_fix_err = 99.0
  b_fix_err = 99.0
  k_fix = 1.0
  k_fix_err = 99.0
  adjust_spectrum,r,b_fix,b_fix_err,r_fix,r_fix_err,k_fix,k_fix_err,stis=stis,quiet=quiet,noplot=noplot
  add_missing_spectrum,r,add=add,compvalues=compvalues,dup=dup,hot=hot,ebv=ebv,noplot=noplot
  if n_elements(compvalues) gt 1 then fitteff = compvalues[6] else fitteff=0d0
  if add eq 1 and quiet eq 0 then print,compvalues[6]
  get_spectrum,r,lambda,spec,err,stis=stis,wein=0

  counter = 0
  while satisfied eq 0 do begin
    ;;offset_thresh = 5 ;; 5
    offset_thresh = 1d30
    if r.name eq 'KIC5608384' or r.name eq 'HD222259B' or r.name eq 'HD222259' or r.name eq 'USco1610' or r.name eq 'LP358-348' or r.name eq 'HD283869' or r.name eq 'EPIC204376071' or r.name eq 'HD221416' or strpos(r.name,'TIC') ne -1 or r.name eq 'HD110082' or nostop eq 1 then offset_thresh = 50000
    counter++
    if counter eq 5 then satisfied = 1
    ;; now everything is in the file, read in that file
    readinphot,r,lambda,spec,err,phot_systems,phot_system,phot_band,phfluxes,phfluxes_err,corr,corr_err,fluxes,fluxes_err,lams,fwhms

    ;;forprint,phfluxes,fluxes

    ;;stop
    ;; calculation of mastercorr and error
    l = where(phot_band ne 'U' and (phot_band ne 'I' or phot_system ne 'Johnson') and phot_system ne 'Wise')
    mastercorr = median(corr[l])
    offset = abs(corr-mastercorr)/sqrt(corr_err^2.0) ;+mastercorr_Err^2.0)
    l = where(phot_band ne 'U' and (phot_band ne 'I' or phot_system ne 'Johnson') and phot_system ne 'Wise' and offset lt offset_thresh,COMPLEMENT=m)
    mastercorr = wmean(corr[l],corr_err[l],error=mastercorr_Err)
    offset = abs(corr-mastercorr)/sqrt(corr_err^2.0) ;+mastercorr_Err^2.0)
    l = where(phot_band ne 'U' and (phot_band ne 'I' or phot_system ne 'Johnson') and phot_system ne 'Wise' and offset lt offset_thresh,COMPLEMENT=m) ;  and phot_system ne 'Hipparcos'
    mastercorr = wmean(corr[l],corr_err[l],error=mastercorr_Err)

    ;; here's where we test to see if there's a problem with blue, red, or NIR regions:
    blue = where(phot_band ne 'U' and phot_band ne 'L' and phot_system ne 'Wise' and offset lt offset_thresh and lams lt 0.52)
    red =  where(phot_band ne 'U'  and (phot_band ne 'I' or phot_system ne 'Johnson') and phot_band ne 'L' and phot_system ne 'Wise' and offset lt offset_thresh and lams gt 0.52 and lams lt 0.9)
    nir =  where(phot_band ne 'L' and phot_system ne 'Wise' and offset lt offset_thresh and lams gt 0.9)
    if nir[0] eq -1 then begin
      print,'failed to find a good NIR match'
      nir =  where(phot_band ne 'L' and phot_system ne 'Wise' and lams gt 0.9)
      stop
    endif
    if blue[0] eq -1 then begin
      mastercorr_b = 1.0
      mastercorr_b_err = 10.
    endif else mastercorr_b = wmean(corr[blue],corr_err[blue],error=mastercorr_b_err)
    if red[0] eq -1 then begin
      mastercorr_r = 1.0
      mastercorr_r_err = 10.
    endif else mastercorr_r = wmean(corr[red],corr_err[red],error=mastercorr_r_err)
    mastercorr_nir = wmean(corr[nir],corr_err[nir],error=mastercorr_nir_err)

    if quiet eq 0 then print,'Blue Correction: '+string(mastercorr_b,format="(D6.3)")+' +/- '+string(mastercorr_b_err,format="(D5.3)")
    if quiet eq 0 then print,'Red Correction:  '+string(mastercorr_r,format="(D6.3)")+' +/- '+string(mastercorr_r_err,format="(D5.3)")
    if quiet eq 0 then print,'NIR Correction:  '+string(mastercorr_nir,format="(D6.3)")+' +/- '+string(mastercorr_nir_err,format="(D5.3)")

    ;; test to see if the errors are inconsistent at 3-sigma:
    r_b = (mastercorr_r-mastercorr_b)/sqrt(mastercorr_b_err^2.0+mastercorr_r_err^2.0)
    nir_r = (mastercorr_nir-mastercorr_r)/sqrt(mastercorr_r_err^2.0+mastercorr_nir_err^2.0)
    if quiet eq 0 then print,'Disagreement: ',r_b,nir_r
    if stis eq 1 then r_b = 0.

    ;if abs(nir_r) gt 3 and r.name ne 'USco1610' and r.name ne 'LP358-348' or r.name eq 'HD283869' then begin
    ;   print,'Fuck'
    ;   stop
    ;endif

    thresh_nir = 4 ;; maximum difference in standard deviations, important parameter
    thresh_b = 3.0 ;; maximum difference in standard deviations, blue
    thresh_nir = 1d5 ;; maximum difference in standard deviations, important parameter
    thresh_b = 1d5 ;; maximum difference in standard deviations, blue
    if r.name eq 'KIC5608384' or r.name eq 'HD222259B' or r.name eq 'HD222259' or r.name eq 'USco1610' or r.name eq 'LP358-348' or r.name eq 'HD283869' or r.name eq 'EPIC204376071' or r.name eq 'HD221416' or strpos(r.name,'TIC') ne -1 or r.name eq 'HD110082' or nostop eq 1 then begin
      thresh_b = 5000
      thresh_nir = 5000
    endif
    if r.name eq 'PM_I16148+6038' then thresh_b = 1.5
    if abs(nir_r) gt thresh_nir or abs(r_b) gt thresh_b then begin
      if abs(nir_r) gt thresh_nir then begin
        r_fix*=(mastercorr_r/mastercorr_nir)
        r_fix_err = r_fix*sqrt((mastercorr_r_err/mastercorr_r)^2.+(mastercorr_nir_err/mastercorr_nir)^2.)
      endif
      if abs(r_b) gt thresh_b then begin
        b_fix*=(mastercorr_b/mastercorr_r)
        b_fix_err = b_fix*sqrt((mastercorr_b_err/mastercorr_b)^2.+(mastercorr_r_err/mastercorr_r)^2.)
      endif
      ll = where(phot_band eq 'K')
      if ll[0] ne -1 then begin
        if n_elements(ll) eq 1 then begin
          k_fix = corr[ll]
          k_fix_err = corr_Err[ll]
        endif else begin
          k_fix = wmean(corr[ll],corr_err[ll],error=k_fix_err)
        endelse
      endif else begin
        k_fix = 1.0
        k_fix_err = 99.0
      endelse
      adjust_spectrum,r,b_fix,b_fix_err,r_fix,r_fix_err,k_fix,k_fix_err,quiet=quiet,noplot=noplot
      add_missing_spectrum,r,add=add,dup=dup,hot=hot,ebv=ebv,noplot=noplot
      get_spectrum,r,lambda,spec,err,stis=stis
    endif else satisfied = 1
    counter++
    if counter gt 10 then begin
      print,'Failed to converge on reasonable solution'
      stop
    endif
    close,/all
  endwhile
  if quiet eq 0 then print,'Final Master correction: '+string(mastercorr,format="(D6.3)")+' +/- '+string(mastercorr_err,format="(D6.3)")

  spec = spec*mastercorr
  err = err*mastercorr

  ;; now fix the RJ tail
  ;; what we want to do is do a RJ fit but also include the photomeric points
  ;; another option is to try extrapolating the models all the way out.
  ;;print,'Fitting Rayleigh-Jeans Tail'
  set_plot,'x'
  !p.multi=[0,1,1]
  q = where(lams gt 2.5)
  if r.name eq 'TYC7577-172-1' then q = where(lams gt 2.0 and lams lt 10)
  if r.name eq 'HD63433' then q = where(lams gt 2.0 and lams lt 30)
  if q[0] ne -1 then begin
    ;; clear out the old shitty RJ fit
    g = where(lambda lt 2 or err ne 0)
    spec = spec[g]
    lambda = lambda[g]
    err = err[g]
    ;; now for a new fit!
    f = where(lambda gt 2.6 and err ne 0 and spec and finite(spec) eq 1)
    if f[0] eq -1 then f = where(lambda gt 2.0 and err ne 0  and finite(spec) eq 1)
    if r.name eq 'TYC7577-172-1' then f = where(lambda gt 3.0 and lambda lt 10)
    if r.name eq 'HD63433' then f = where(lambda gt 7.0 and lambda lt 10)
    guess = [1d-10,5]
    result = mpfitfun('MYRJL',[lambda[f],lams[q]],[spec[f],phfluxes[q]],[spec[f]/2.,phfluxes_err[q]],guess,perror=perr,/quiet)
    newlambda = generatearray(10.0,30.0,100)
    newspec = myrjl(newlambda,result)
    newerr = 0.0*((newspec)*(perr[0]/result[0]))

    lambda = [lambda,newlambda]
    spec = [spec,newspec]
    err = [err,newerr]
    if noplot eq 0 then begin
      plot,lambda,spec,xrange=[2,30],/xlog,/ylog,/xstyle
      oploterror,lams[q],phfluxes[q],fwhms[q],phfluxes_err[q],psym=8,color=cgcolor('red')
    endif

    off = dblarr(n_elements(q))
    for jj = 0,n_elements(q)-1 do begin
      get_flux,phot_systems,phot_system[q[jj]],phot_band[q[jj]],lambda,spec,err,flux,flux_err

      if noplot eq 0 then oplot,[lams[q[jj]]],[flux],psym=8
      off[jj] = (phfluxes[q[jj]]-flux)/phfluxes_err[q[jj]]
    endfor
    f = where(lambda gt 2.6 and err ne 0)
    guess = [1.0,-1.0]
    newlambda = generatearray(max(lambda[where(err ne 0)]),30.0,100)
    nmonte = 500
    bigarr = dblarr(n_elements(newlambda),nmonte)
    newerr = dblarr(n_elements(newlambda))
    for iii = 0,nmonte-1 do begin
      res_temp = result+perr*randomn(seed,2)
      bigarr[*,iii] = MYRJL(newlambda,res_temp)
    endfor
    for jjj = 0,n_elements(newlambda)-1 do begin
      newerr[jjj] = stdev(bigarr[jjj,*])
    endfor
    g = where(err eq 0 and lambda gt 2.0)
    while n_elements(newerr) lt n_elements(g) do newerr = [newerr,median(newerr)]
    while n_elements(newerr) gt n_elements(g) do g = [min(g)-1,g]
    err[g] = newerr
  endif

  ;; re-extract all synthetic photometry:
  fluxes = dblarr(n_elements(phot_system))
  fluxes_err = fluxes
  for jj = 0,n_elements(phot_system)-1 do begin
    get_flux,phot_systems,phot_system[jj],phot_band[jj],lambda,spec,err,flux,flux_err
    fluxes[jj] = flux
    fluxes_err[jj] = flux_err
  endfor
  readinphot,r,lambda,spec,err,phot_systems,phot_system,phot_band,phfluxes,phfluxes_err,corr,corr_err,fluxes,fluxes_err,lams,fwhms
  ;;qq = where(finite(corr_err) eq 0)
  ;;if qq[0] ne -1 then corr_err[qq] = corr[qq]

  ;; calculation of mastercorr and error
  test_corr = median(corr)
  offset = abs(test_corr-corr)/sqrt(corr_err^2.0) ;+mastercorr_Err^2.0)
  l = where(offset lt offset_thresh and corr_err lt 1,COMPLEMENT=m)
  test_corr2 = wmean(corr[l],corr_err[l],error=test_corr_Err)
  offset = abs(test_corr2-corr)/sqrt(corr_err^2.0) ;+mastercorr_Err^2.0)
  l = where(offset lt offset_thresh,COMPLEMENT=m)
  l2 = where(offset lt offset_thresh,COMPLEMENT=m)
  l_tmp = where((offset lt offset_thresh and offset gt 0) or phot_system eq 'Wise',COMPLEMENT=m2)
  l_nowise = where(offset lt offset_thresh and phot_system ne 'Wise')
  l_cover = where(offset lt offset_thresh and phot_system ne 'Wise' and lams-fwhms gt 0.31)
  if wexcess eq 1 then begin
    l = where(offset lt offset_thresh and lams lt 10,COMPLEMENT=m)
    l2 = where(offset lt offset_thresh,COMPLEMENT=m)
  endif

  rchisq = total((corr[l]-test_corr2)^2.0/corr_err[l]^2.0)/(n_elements(l)-2)

  rchisq_nowise = total((corr[l_nowise]-test_corr2)^2.0/corr_err[l_nowise]^2.0)/(n_elements(l_nowise)-2)
  rchisq_cover = total((corr[l_cover]-test_corr2)^2.0/corr_err[l_cover]^2.0)/(n_elements(l_cover)-2)
  forprint, string(phot_system[l2],format="(A10)")+string(9b)+string(phot_band[l2],format="(A4)")+string(9b)+string(lams[l2],format="(D5.2)"),$
  (corr[l2]-test_corr)/corr_err[l2], corr[l2],corr_err[l2],textout='phots/'+name+'.phot',comment='-----------------Goddies-----------------',/silent
  if m[0] ne -1 then begin
    forprint,string(phot_system[m],format="(A10)")+string(9b)+string(phot_band[m],format="(A4)")+string(9b)+string(lams[m],format="(D5.2)"),$
    (corr[m]-test_corr)/corr_err[m],(corr[m]),corr_err[m],textout='phots/'+name+'.badphot',comment='-----------------Baddies-----------------',/silent
    if quiet eq 0 then begin
      print,phot_system[m]
      print,phot_band[m]
    endif
    ;;stop
  endif
  if quiet eq 0 then print,'Final Secondary Correction: '+string(test_corr2,format="(D6.3)") +' +/- '+string(test_corr_err,format="(D6.3)")
  if m[0] eq -1 then count = 0 else count = n_Elements(m)
  if m2[0] eq -1 then count2 = 0 else count2 = n_Elements(m2)
  if quiet eq 0 then print,'Final Reduced Chi^2: '+string(rchisq,format="(D5.2)")+' with '+string(n_elements(l),format="(I3)")+' photometric points.  '+string(count,format="(I3)")+' points excluded. '+string(count2,format="(I2)")+' non-Wise non-I,R points excluded.'
  if quiet eq 0 then print,'Final Reduced Chi^2 without wise: '+string(rchisq_nowise,format="(D5.2)")
  if quiet eq 0 then print,'Final Reduced Chi^2 points covered: '+string(rchisq_cover,format="(D5.2)")
  goodpoints = n_elements(l_nowise)

  ;;;
  if test_corr2 gt 0 then spec*=test_corr2
  if r.cns3 eq 'GJ 380' or r.cns3 eq 'GJ 436' or r.cns3 eq 'GJ 725A' then spec*=1.04 ;; THESE NEED FIXING
  if r.cns3 eq 'GJ 380' then spec*=1.01

  gg = where(err ne 0)
  medsnr = median(spec[gg]/err[gg])
  loc = where(finite(spec) eq 1 and finite(lambda) eq 1)
  fbol = integral(lambda[loc],spec[loc])*1d12

  nmonte = 200
  fbols = dblarr(nmonte)
  for kk = 0,nmonte-1 do begin
    tmp_spec = spec+err*randomn(seed,n_elements(spec))
    tmp_spec*=(1.+(randomn(seed)*test_corr_err))
    jj = where(finite(tmp_spec) eq 1 and finite(lambda) eq 1)
    fbols[kk] = integral(lambda[jj],tmp_spec[jj])
  endfor
  fbol_err = stdev(fbols)*1d12
  ;; add 0.5% systematic error
  fbol_err=sqrt(fbol_err^2.+(0.005*fbol)^2.)
  if quiet eq 0 then print,'Final Fbol: '+string(fbol,format="(D9.4)")+' +/- '+string(fbol_err,format="(D9.4)")
  close,/all


  ;; Plot!
  ;;if add eq 1 then begin
  gcal_plotter2,r,phot_system,offset,offset_thresh,name,lambda,spec,err,lams,phfluxes,phfluxes_err,fwhms,$
  FLUXES,FLUXES_err,add=add,pxrange=pxrange,stis=stis,photometrycolor=photometrycolor,nolegend=nolegend
  ;;gcal_plotter,r,phot_system,offset,offset_thresh,name,lambda,spec,err,lams,phfluxes,phfluxes_err,fwhms,FLUXES,FLUXES_err,add=add
  close,/all
  ;;endif

  ;; save the adjusted spectrum
  finalspec = spec
  finallambda = lambda
  finalerr = err
  finalfbol = fbol
  finalfbol_err = fbol_err


END



PRO adjust_spectrum,h,op_adjust,op_adjust_err,IR_adjust,ir_adjust_err,k_fix,k_fix_err,stis=stis,quiet=quiet,noplot=noplot

  if n_elements(k_fix) eq 0 or n_elements(k_fix_err) eq 0 then begin
    k_fix = 1.0
    k_fix_err = 99.
  endif
  COMPILE_OPT idl2, HIDDEN
  set_plot,'x'
  !p.multi=[0,1,2]
  snifsfwhm = (7000.0/1000.0)
  c = 299792.458D

  extractdata,h,opspec,operr,oplambda,opspec_b,operr_b,oplambda_b,oprv,irspec,irerr,irlambda,irrv,stis=stis

  if h.name eq 'GJ673' then irspec[where(irlambda gt 1.9)]*=0.93
  if k_fix ne 1. and abs((k_fix-1.0)/k_fix_err) gt 1 then begin
    if k_fix gt 1.02 then k_fix = 1.02
    if k_fix lt 0.98 then k_fix = 0.98
    irspec[where(irlambda gt 1.9)]*=k_fix
  endif

  ;; apply RV corrections
  irlambda = irlambda/(irrv/c+1.)

  opoverlaplambda = oplambda[where(oplambda gt min(irlambda))]
  iroverlaplambda = irlambda[where(irlambda lt max(oplambda))]
  opoverlap = opspec[where(oplambda gt min(irlambda))]
  iroverlap = irspec[where(irlambda lt max(oplambda))]
  opoverlap_err = operr[where(oplambda gt min(irlambda))]
  iroverlaplambda_err = irerr[where(irlambda lt max(oplambda))]

  compspec = interpol(iroverlap,iroverlaplambda,opoverlaplambda) ;; for now we interpolate the IR to the optical
  l = opoverlaplambda
  limits = [0.805,0.90]
  rv = 0d0
  if h.name eq 'GJ_105A' then limits = [0.92,0.94]
  if n_elementS(opoverlap) gt 5 then compute_xcor_rv, l, opoverlap,'M1.5', rv, wavelim=[min(l),max(l)],norm_reg=limits,twave=l, tflux=compspec, maxshift=20,showplot=0
  if h.name eq 'LSPM_J0315+0103' then rv = 0.0
  if h.name eq 'NLTT37349' then rv = 0.0
  oplambda = oplambda/(rv/c+1.)
  if stis eq 0 then oplambda_b = oplambda_b/(rv/c+1.)



  limits = [0.805,0.88]


  right = irspec[where(irlambda gt limits[0] and irlambda lt limits[1])]
  right_err = irerr[where(irlambda gt limits[0] and irlambda lt limits[1])]
  right_lambda = irlambda[where(irlambda gt limits[0] and irlambda lt limits[1])]
  if n_elements(right) eq 1 then s1 = right else s1 = median(right)
  left = opspec[where(oplambda gt limits[0] and oplambda lt limits[1])]
  left_err = operr[where(oplambda gt limits[0] and oplambda lt limits[1])]
  left_lambda = oplambda[where(oplambda gt limits[0] and oplambda lt limits[1])]
  right = interpol(right,right_lambda,left_lambda)
  s2 = median(left)
  change = right/left
  coverlap_base = median(change)
  coverlap_base_err = sqrt(stdev(change)^2.0+(0.02*coverlap_base)^2.0) ;; 1-2% is the error in the flux cal
  coverlap_other = coverlap_base*ir_adjust
  coverlap_other_err = ir_adjust_err*coverlap_base
  if ir_adjust_err lt 1 then coverlap = wmean([coverlap_base,coverlap_other],[coverlap_base_err,coverlap_other_err]) else coverlap = coverlap_base

  if finite(coverlap) eq 0 then coverlap = 1
  opspec*=Coverlap
  operr*=Coverlap
  if quiet eq 0 then print,'Red correction  ',coverlap

  if stis eq 0 then begin
    ;; same for blue/red channel
    limits2 = [0.510,0.515]
    if h.name eq 'PM_I13457+1453' then limits = [0.50,0.511]
    if h.name eq 'PM_I04073-2429' then limits2 = [0.5,0.512]
    right = opspec[where(oplambda gt limits2[0] and oplambda lt limits2[1])]
    right_err = operr[where(oplambda gt limits2[0] and oplambda lt limits2[1])]
    right_lambda = oplambda[where(oplambda gt limits2[0] and oplambda lt limits2[1])]
    s1 = median(right)
    left = opspec_b[where(oplambda_b gt limits2[0] and oplambda_b lt limits2[1])]
    left_err = operr_b[where(oplambda_b gt limits2[0] and oplambda_b lt limits2[1])]
    left_lambda = oplambda_b[where(oplambda_b gt limits2[0] and oplambda_b lt limits2[1])]
    right = interpol(right,right_lambda,left_lambda)
    s2 = median(left)
    Coverlap = s1/s2

    change = right/left
    coverlap_base = median(change)
    coverlap_base_err = sqrt(stdev(change)^2.0+(0.03*coverlap_base)^2.0) ;; 1-2% is the error in the flux cal
    coverlap_other = coverlap_base*op_adjust
    coverlap_other_err = op_adjust_err*coverlap_base
    if op_adjust_err lt 1 then coverlap = wmean([coverlap_base,coverlap_other],[coverlap_base_err,coverlap_other_err]) else coverlap = coverlap_base

    ;;coverlap*=op_adjust
    if stis eq 0 then begin
      opspec_b*=Coverlap
      operr_b*=Coverlap
    endif
  endif
  if noplot eq 0 then begin
    plot,oplambda,opspec,xrange=[0.475,0.55]
    if stis eq 0 then begin
      oplot,oplambda_b,opspec_b,color=cgcolor('blue')
      if quiet eq 0 then print,'Blue correction  ',coverlap
    endif
  endif

  limits = [0.83,0.9]                              ;;[0.925,0.94];;

  opspec_o = opspec[where(oplambda gt limits[0] and oplambda lt limits[1])]
  operr_o = operr[where(oplambda gt limits[0] and oplambda lt limits[1])]
  oplambda_o = oplambda[where(oplambda gt limits[0] and oplambda lt limits[1])]
  irspec_o = irspec[where(irlambda gt limits[0] and irlambda lt limits[1])]
  irlambda_o = irlambda[where(irlambda gt limits[0] and irlambda lt limits[1])]
  irerr_o = irerr[where(irlambda gt limits[0] and irlambda lt limits[1])]

  ;;irspec_o = gaussfold(irlambda_o,irspec_o,snifsfwhm) ; convolve to expected resolution
  irspec_o = interpol(irspec_o,irlambda_o,oplambda_o)
  numbefore = 1.0*n_elements(irerr_o)
  irerr_o = interpol(irerr_o,irlambda_o,oplambda_o)
  numbafter = 1.0*n_elements(irerr_o)
  irerr_o*=sqrt(numbafter/numbefore) ; we gain Signal by root(N) when we bin up data
  toterr = (irerr_o + operr_o)
  overlap = (irspec_o*(operr_o/toterr) + opspec_o*(irerr_o/toterr))
  overlap_err = 1.0/sqrt(1.0/(irerr_o^2.0) + 1.0/(operr_o^2.0))
  if h.name eq 'PM_I11033+3558' then begin
    overlap = irspec_o
    overlap_err = irerr_o
  endif

  if stis eq 0 then begin
    masterspec = [opspec_b[where(oplambda_b lt min(oplambda))],opspec[where(oplambda lt limits[0])],overlap,irspec[where(irlambda gt limits[1])]]
    masterlambda = [oplambda_b[where(oplambda_b lt min(oplambda))],oplambda[where(oplambda lt limits[0])],oplambda_o,irlambda[where(irlambda gt limits[1])]]
    mastererr = [operr_b[where(oplambda_b lt min(oplambda))],operr[where(oplambda lt limits[0])],overlap_err,irerr[where(irlambda gt limits[1])]]
  endif else begin
    masterspec = [opspec[where(oplambda lt limits[0])],overlap,irspec[where(irlambda gt limits[1])]]
    masterlambda = [oplambda[where(oplambda lt limits[0])],oplambda_o,irlambda[where(irlambda gt limits[1])]]
    mastererr = [operr[where(oplambda lt limits[0])],overlap_err,irerr[where(irlambda gt limits[1])]]
  endelse

  l = where(finite(masterspec) eq 1 and finite(masterlambda) eq 1 and finite(mastererr) eq 1 and mastererr gt 0)
  masterspec = masterspec[l]
  masterlambda = masterlambda[l]
  mastererr = mastererr[l]
  if noplot eq 0 then begin
    plot,masterlambda,masterspec,/xstyle,/ystyle,xrange=[0.68,1.1],thick=1
    oplot,oplambda,opspec,color=cgcolor('green')
    oplot,irlambda,irspec,color=cgcolor('red')
  endif


  h.fullspec = ptr_new(masterspec)
  h.fulllambda = ptr_new(masterlambda)
  h.fullerr = ptr_new(mastererr)


END


PRO extractdata,h,opspec,operr,oplambda,opspec_b,operr_b,oplambda_b,oprv,irspec,irerr,irlambda,irrv,stis=stis

  stis = 0
  COMPILE_OPT idl2, HIDDEN

  if h.name eq 'HD29391' then begin
    opspec = *h.FULLSPEC
    oplambda = *h.FULLLAMBDA*10000d0
    operr = *h.FULLERR
  endif else begin
    if h.uh22snr gt 50 then begin
      operr = 1
      get_struct_wave_flux, h, oplambda, opspec, var=operr, /restwave
      operr = sqrt(operr)
      operr_b = 1
      get_struct_wave_flux, h, oplambda_b, opspec_b, var=operr_b, /restwave, /blue
      operr_b = sqrt(operr_b)
    endif else begin ;; get a template optical spectrum
      stop
      spawn,'ls pickles/*.dat',list
      spt = strlowcase(h.spt)
      qq = where(strpos(list,spt) ne -1)
      if qq[0] eq -1 then stop
      readcol,list[qq[0]],oplambda,opspec
      gg = wherE(oplambda lt 1d4)
      oplambda = oplambda[gg]
      opspec = opspec[gg]
      operr = opspec*0.02
      stis = 1
    endelse
  endelse


  if h.stisspec ne ptr_new() then begin
    tmp = *h.stisspec
    if n_elements(tmp) gt 10 then begin
      oplambda = tmp[*,0]
       opspec = tmp[*,1]
       operr = tmp[*,2]
       stis = 1
    endif
  endif

  oplambda /= 10000.0
  if stis ne 1 then oplambda_b /= 10000.0

  if h.irtfsnr gt 10 then begin
    temp = *h.irtfspec
    if (size(temp))[0] eq 2 then begin
      irlambda = temp[*,0]
      irspec = temp[*,1]
      irerr = temp[*,2]
      irrv = h.irrv
    endif else begin
      irlambda = [temp[*,0,5],temp[*,0,4],temp[*,0,3],temp[*,0,2],temp[*,0,1],temp[*,0,0]]
      irspec = [temp[*,1,5],temp[*,1,4],temp[*,1,3],temp[*,1,2],temp[*,1,1],temp[*,1,0]]
      irerr = [temp[*,2,5],temp[*,2,4],temp[*,2,3],temp[*,2,2],temp[*,2,1],temp[*,2,0]]
      irrv = h.irrv
    endelse
  endif else begin ;; get a template optical spectrum
    spawn,'ls pickles/*.dat',list
    spt = strlowcase(h.spt)
    qq = where(strpos(list,spt) ne -1)
    if qq[0] eq -1 then stop
    irlambda = 0
    f = 0
    while max(irlambda) lt 2.0 do begin
      readcol,list[qq[f]],irlambda,irspec
      gg = wherE(irlambda gt 0.8)
      irlambda = irlambda[gg]
      irspec = irspec[gg]
      irlambda/=1d4
      irerr = irspec*0.02
      irrv = 0
      f++
    endwhile
  endelse

  q = where(irlambda lt 2.54)
  irlambda = irlambda[q]
  irspec = irspec[q]
  irerr = irerr[q]


  order = sort(irlambda)
  irlambda = irlambda[order]
  irspec = irspec[order]
  irerr = irerr[order]
  locs = where(finite(irspec) eq 1 and finite(irerr) eq 1)
  irlambda = irlambda[locs]
  irspec = irspec[locs]
  irerr = irerr[locs]

END



FUNCTION myRJL,X,P
  COMPILE_OPT idl2, HIDDEN
  t = p[0]/x^p[1]               ; + p[1]/x^4.0
  RETURN, t
END

FUNCTION myWein,X,P
  COMPILE_OPT idl2, HIDDEN
  t = (p[0]/x^5.0)*exp(p[1]/x)
  RETURN, t
END


;; used to combine 2 spectra
PRO combine,c1,c2,arr

  tmp1 = *c1.irtfspec
  l1 = tmp1[*,0]
  s1 = tmp1[*,1]
  err1 = tmp1[*,2]
  tmp2 = *c2.irtfspec
  l2 = tmp2[*,0]
  s2 = tmp2[*,1]
  err2 = tmp2[*,2]
  s2 = interpol(s2,l2,l1)
  s = (s1/median(s1)+s2/median(s2))/2.
  arr = [[l1],[s],[err1]]

END

