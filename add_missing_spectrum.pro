PRO add_missing_spectrum,r,add=add,compvalues=compvalues,dup=dup,hot=hot,conroy=conroy,plot=plot,modelblue=modelblue,ebv=ebv,noplot=noplot

  if n_elements(dup) eq 0 then dup = 0
  if n_elements(ebv) eq 0 then ebv = 0
  if n_elements(modelblue) eq 0 then modelblue = 0
  if n_elements(plot) eq 0 then plot = 1
  if n_elements(conroy) eq 0 then conroy = 0
  if n_elements(add) eq 0 then add = 1
  if n_elements(hot) eq 0 then hot = 0
  if n_elements(noplot) eq 0 then noplot = 0
  if noplot eq 1 then plot = 0
  if add ge 1 then begin
     COMPILE_OPT IDL2, HIDDEN
     models = mrdfits('Models_CIFIST_Aug2013_filler.fits',1,/silent)

     ;; if ebv > 0 then redden the models
     modelspectra = models.spectrum
     modelheader = models.header
     lambda_m = sxpar(modelheader,'LAMBDA_0') + findgen(sxpar(modelheader,'NLAMBDA'))*sxpar(modelheader,'D_LAMBDA')
     teff = models.teff
     logg = abs(models.logg)
     metal = models.metal
     if ebv>0 then begin
         redmodelspectra = modelspectra*0
         for i=0,n_elements(modelspectra[*,0])-1 do begin
            ccm_unred,lambda_m,modelspectra[i,*],-1d0*ebv,redspec
            redmodelspectra[i,*] = redspec
            ; plot,lambda_m,modelspectra[i,*],/xstyle,/ystyle,/xlog,/ylog
            ; oplot,lambda_m,redspec
            ; stop
         endfor
         models.spectrum = redmodelspectra
         modelspectra = redmodelspectra
     endif
     
     mh = r.irtffeh
     if mh lt -1.0 then metalrange=[-1.6,-0.9]
     if mh ge -1.0 and mh lt -0.5 then metalrange=[-1.1,-0.4]
     if mh ge -0.5 and mh lt 0.0 then metalrange=[-0.6,0.1]
     if mh ge 0.0 and mh lt 0.3 then metalrange=[-0.1,0.4]
     if mh ge 0.3 then metalrange=[0.2,0.6]
     if n_elements(dup) gt 0 then metalrange = [-0.3,0.3]
     
     flux = *r.fullspec
     wave = *r.fulllambda*1d4
     var = (*r.fullerr)^2.0

     max_lambda = 1d5
     l = where(wave gt 10000)
     fwhm = 7000./1300.
     if n_elements(dup) gt 0 then fwhm = median(wave-shift(wave,1))*2.0
     run_j9_on_models_2,wave[l],flux[l],var[l],models,compvalues, metalrange=metalrange, savespec=savespec,savemetal=savemetal,saveteff=saveteff,savelogg=savelogg,saveweight=saveweight,dup=dup,fwhm=fwhm,plot=0,hot=hot
          
     flux = *r.fullspec
     model1 = modelspectra[(where(teff eq saveteff[0] and abs(logg) eq abs(savelogg[0]) and metal eq savemetal[0]))[0],*]
     model2 = modelspectra[(where(teff eq saveteff[1] and abs(logg) eq abs(savelogg[1]) and metal eq savemetal[1]))[0],*]
     model3 = modelspectra[(where(teff eq saveteff[2] and abs(logg) eq abs(savelogg[2]) and metal eq savemetal[2]))[0],*]
     mspec = saveweight[0]*model1/mean(model1[where(finite(model1) eq 1)]) + (saveweight[1]-saveweight[0])*model2/mean(model2[where(finite(model2) eq 1)]) + (1.-saveweight[1])*model3/mean(model3[where(finite(model3) eq 1)])
     if n_elements(dup) gt 0 then mspec = gaussfold(lambda_m,mspec,fwhm) else  mspec = gaussfold(lambda_m,mspec,fwhm)
     
     modelspec = mspec*(median(flux[where(wave gt 1.7d4 and wave lt 1.75d4)])/median(mspec[where(lambda_m gt 1.7d4 and lambda_m lt 1.75d4)]))
     dlambda = abs(median(wave-shift(wave,1)))
     if n_elements(dup) gt 0 then dlambda = 20.
     if r.name eq 'SD0423-04' or r.name eq '2M1728+39' or r.name eq '2M2252-17' then dlambda = 50.
     intermodelwave = generatearray(500,max_lambda-50,(max_lambda)/dlambda)
     modelspec_base = modelspec
     modelspec = interpol(modelspec,lambda_m,intermodelwave)
     w = wave/1d4

     ;;set_plot,'x'
     ;;print,max_lambda
     ;;plot,intermodelwave,modelspec,/xlog,/ylog,/xstyle,/ystyle
     ;;plot,lambda_m,modelspec_base,/xlog,/ylog,/xstyle,/ystyle

     if conroy eq 0 then begin
        pl1 = 1.34
        pl2 = 1.45
        pl3 = 1.77
        pl4 = 2.1
     endif else begin
        pl1 = 1.35
        pl2 = 1.43
        pl3 = 1.79
        pl4 = 2.1
     endelse

     
     if add eq 1 then begin
        if modelblue eq 1 then begin
           fancyflux = [modelspec[where(intermodelwave lt min(wave))],$
                        flux[where(w le pl1)],modelspec[where(intermodelwave gt pl1*1d4 and intermodelwave le pl2*1d4)],$
                        flux[where(w gt pl2 and w le pl3)],modelspec[where(intermodelwave gt pl3*1d4 and intermodelwave le pl4*1d4)],$
                        flux[where(w gt pl4)],modelspec[where(intermodelwave gt max(wave))]]
           fancywave = [intermodelwave[where(intermodelwave lt min(wave))],$
                        wave[where(w le pl1)],intermodelwave[where(intermodelwave gt pl1*1d4 and intermodelwave le pl2*1d4)],$
                        wave[where(w gt pl2 and w le pl3)],intermodelwave[where(intermodelwave gt pl3*1d4 and intermodelwave le pl4*1d4)],$
                        wave[where(w gt 2.1)],intermodelwave[where(intermodelwave gt max(wave))]]
           fancyerr = [dblarr(n_elements(where(intermodelwave lt min(wave)))),$
                       var[where(w le pl1)],dblarr(n_elements(where(intermodelwave gt pl1*1d4 and intermodelwave le pl2*1d4))),$
                       var[where(w gt pl2 and w le pl3)],dblarr(n_elements(where(intermodelwave gt pl3*1d4 and intermodelwave le pl4*1d4))),$
                       var[where(w gt pl4)],dblarr(n_elements(where(intermodelwave gt max(wave))))]
        endif else begin
           fancyflux = [$       ;modelspec[where(intermodelwave lt min(wave))],$
                       flux[where(w le pl1)],modelspec[where(intermodelwave gt pl1*1d4 and intermodelwave le pl2*1d4)],$
                       flux[where(w gt pl2 and w le pl3)],modelspec[where(intermodelwave gt pl3*1d4 and intermodelwave le pl4*1d4)],$
                       flux[where(w gt pl4)],modelspec[where(intermodelwave gt max(wave))]]
           fancywave = [$       ;intermodelwave[where(intermodelwave lt min(wave))],$
                       wave[where(w le pl1)],intermodelwave[where(intermodelwave gt pl1*1d4 and intermodelwave le pl2*1d4)],$
                       wave[where(w gt pl2 and w le pl3)],intermodelwave[where(intermodelwave gt pl3*1d4 and intermodelwave le pl4*1d4)],$
                       wave[where(w gt 2.1)],intermodelwave[where(intermodelwave gt max(wave))]]
           fancyerr = [$        ;dblarr(n_elements(where(intermodelwave lt min(wave)))),$
                      var[where(w le pl1)],dblarr(n_elements(where(intermodelwave gt pl1*1d4 and intermodelwave le pl2*1d4))),$
                      var[where(w gt pl2 and w le pl3)],dblarr(n_elements(where(intermodelwave gt pl3*1d4 and intermodelwave le pl4*1d4))),$
                      var[where(w gt pl4)],dblarr(n_elements(where(intermodelwave gt max(wave))))]
        endelse
        if n_elements(dup) gt 0 then begin
           redcut3 = 2.0
           bluecut3 = 1.8
           llm = where(intermodelwave lt min(Wave)+100d0 and intermodelwave gt min(wave))
           lls = where(wave lt min(wave)+100d0)
           bcoeff = median(flux[lls])/median(modelspec[llm])
           llm = where(intermodelwave gt max(Wave)-100d0 and intermodelwave lt max(wave))
           lls = where(wave gt max(wave)-100d0)
           ;if r.name eq '2M0700+31' then begin
           ;   tmp = findgen(n_elements(where(intermodelwave gt 3d4)))
           ;   modelspec[where(intermodelwave gt 3d4)]*=(tmp/median(tMp))^0.01
           ;endif
           rcoeff = median(flux[lls])/median(modelspec[llm])
           fancyflux = [modelspec[where(intermodelwave lt min(wave))]*bcoeff,$
                       flux[where(w le pl1)],modelspec[where(intermodelwave gt pl1*1d4 and intermodelwave le pl2*1d4)],$
                       flux[where(w gt pl2 and w le bluecut3)],modelspec[where(intermodelwave gt bluecut3*1d4 and intermodelwave le redcut3*1d4)],$
                       flux[where(w gt redcut3)],modelspec[where(intermodelwave gt max(wave))]*rcoeff]
           fancywave = [intermodelwave[where(intermodelwave lt min(wave))],$
                       wave[where(w le pl1)],intermodelwave[where(intermodelwave gt pl1*1d4 and intermodelwave le pl2*1d4)],$
                       wave[where(w gt pl2 and w le bluecut3)],intermodelwave[where(intermodelwave gt bluecut3*1d4 and intermodelwave le redcut3*1d4)],$
                       wave[where(w gt redcut3)],intermodelwave[where(intermodelwave gt max(wave))]]
           fancyerr = [dblarr(n_elements(where(intermodelwave lt min(wave)))),$
                      var[where(w le pl1)],dblarr(n_elements(where(intermodelwave gt pl1*1d4 and intermodelwave le pl2*1d4))),$
                      var[where(w gt pl2 and w le bluecut3)],dblarr(n_elements(where(intermodelwave gt bluecut3*1d4 and intermodelwave le redcut3*1d4))),$
                      var[where(w gt redcut3)],dblarr(n_elements(where(intermodelwave gt max(wave))))]
        endif
           
     endif

     if plot eq 1 then begin
        set_plot,'x'
        !p.multi = [0,1,1,0,0]
     endif
     if r.name eq 'EPIC_205117205' or r.name eq 'EPIC_205046529' then begin
        models2 = mrdfits('~/Dropbox/Radii/Models_CIFIST_Apr2016_long.fits',1,/silent)
        modelspectra2 = models2.spectrum
        modelheader2 = models2.header
        lambda_m2 = sxpar(modelheader2,'LAMBDA_0') + findgen(sxpar(modelheader2,'NLAMBDA'))*sxpar(modelheader2,'D_LAMBDA')
        teff2 = models2.teff
        logg2 = abs(models2.logg)
        metal2 = models2.metal
        newmodel = modelspectra2[where(teff2 eq 3500),*]
        newmodel = newmodel[*]
        qq1 = where(fancywave gt 9d4 and fancywave lt 10d4)
        qq2 = where(lambda_m2 gt 9d4 and lambda_m2 lt 10d4)
        newmodel/=median(newmodel[qq2])/median(fancyflux[qq1])
        
        if plot eq 1 then begin
           plot,w,flux,xrange=[2,30],/xlog,/ylog,yrange=[1d-20,1d-14],/xstyle
           oplot,intermodelwave/1d4,modelspec,color=cgcolor('red')
           past = where(lambda_m2 gt max(fancywave))
           oplot,lambda_m2[past]/1d4,newmodel[past],color=cgcolor('teal')
        endif
        
        fancyflux = [fancyflux,newmodel[past]]
        fancywave = [fancywave,lambda_m2[past]]
        fancyerr = [fancyerr,dblarr(n_elements(past))]
     endif
     
     
     if add eq 2 then begin
        fancyflux = [flux,modelspec[where(intermodelwave gt max(wave))]]
        fancywave = [wave,intermodelwave[where(intermodelwave gt max(wave))]]
        fancyerr = [var,dblarr(n_elements(where(intermodelwave gt max(wave))))]
     endif
     if plot eq 1 then begin
        !p.multi=[0,2,2,0,0]
        plot,fancywave,fancyflux,/xstyle,/ystyle,yrange=[min(fancyflux[where(fancyflux gt 0)]),max(fancyflux)],/xlog,/ylog
        oplot,fancywave[where(fancyerr eq 0.0)],fancyflux[where(fancyerr eq 0.0)],color=cgcolor('red')

        ll = where(fancyflux gt 1d-20)
        plot,fancywave[ll],fancyflux[ll],/xstyle,/ystyle,xrange=[1.2d4,1.55d4],/ylog
        oplot,fancywave[where(fancyerr eq 0.0)],fancyflux[where(fancyerr eq 0.0)],color=cgcolor('red')

        ;;plot,fancywave,fancyflux,/xstyle,/ystyle,xrange=[1.3d4,1.6d4]
        ;;oplot,fancywave[where(fancyerr eq 0.0)],fancyflux[where(fancyerr eq 0.0)],color=cgcolor('red')

        plot,fancywave,fancyflux,/xstyle,/ystyle,xrange=[1.6d4,2.2d4]
        oplot,fancywave[where(fancyerr eq 0.0)],fancyflux[where(fancyerr eq 0.0)],color=cgcolor('red')

        plot,fancywave,fancyflux,/xstyle,/ystyle,xrange=[2.0d4,2.7d4]
        oplot,fancywave[where(fancyerr eq 0.0)],fancyflux[where(fancyerr eq 0.0)],color=cgcolor('red')
        ;;print,r.newteff,compvalues[6],compvalues[7]

        ll1 = where(fancyerr eq 0.0 and fancywave gt 1.3d4 and fancywave lt 1.5d4)
        ll2 = where(fancyerr eq 0.0 and fancywave gt 1.7d4 and fancywave lt 2.2d4)
     endif
     fancydata = [[fancywave],[fancyflux],[fancyerr]]
     ;;stop

  endif else begin
     flux = *r.fullspec
     wave = *r.fulllambda*1d4
     var = (*r.fullerr)^2.0
     fancydata = [[wave],[flux],[var]]
  endelse
  r.fspec = ptr_new(fancydata)

END

