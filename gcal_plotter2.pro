PRO gcal_plotter2,r,phot_system,offset,offset_thresh,name,lambda,spec,sp_err,lams,phfluxes,phfluxes_err,$
                  fwhms,FLUXES,FLUXES_err,photometrycolor=photometrycolor,add=add,inset=inset,long=long,pxrange=pxrange,stis=stis,soar=soar,$
                  oc=oc,e_corr=e_corr,corr=corr,dup=dup,printname=printname,$
                  nolegend=nolegend ;; 0 (default) o-c in terms of sigmas, 1 = in terms of magnitudes

  oc = 0
  fac = 10000d0
  spec*=lambda*fac
  ;;sp_err*=lambda*fac
  phfluxes*=lams*fac
  phfluxes_err*=lams*fac
  FLUXES*=lams*fac
  FLUXES_err*=lams*fac

  soar = 1
  err = sp_err
  if n_elements(nolegend) eq 0 then nolegend = 0
  if n_elements(photometrycolor) eq 0 then photometrycolor = 0
  if n_elements(soar) eq 0 then soar = 0
  if n_elements(oc) eq 0 then oc = 0
  if n_elements(stis) eq 0 then stis = 0
  if n_elementS(printname) eq 0 then printname = 0
  if n_elements(long) eq 0 then long = 0
  ;;if r.name eq 'USco1610' then long = 1
  if n_elements(add) eq 0 then add = 0
  if n_elements(inset) eq 0 then inset = 0
  charsize = 1.25
  thick = 4.0
  xthick = 7
  ythick = 7
  charthick = 6.0
  errthick = 6.0
  plotsym,0,/fill
  if n_elements(pxrange) ne 2 then begin
     PXRANGE = [0.4,10.0]
     xrange = [0.4,10.0]
     xrange_blue = [0.25,0.6]
     xrange_red = [2.0,30]
     if long eq -1 then xrange = [0.45,2.4]
  endif else xrange = pxrange
  yrange = [min(spec[where(lambda gt xrange[0] and lambda lt xrange[1] and spec gt 0)]),max(spec)*1.3]
  if yrange[0] le 0 then yrange[0] = min(spec[wherE(spec gt 1d-20)])
  if yrange[0] le 0 then stop
  if inset eq 1 then thick = 4
  !p.font=0

  angstrom = cgsymbol("Angstrom")
  angstrom ='A'
  sigma = greek('sigma')+'!3'   ;textoidl('\sigma')

  errsym = 8
  errthick=xthick
  l = where(offset lt offset_thresh,COMPLEMENT=m)
  factor = 1d0
  set_plot,'PS'
  ;; colors
  photcolor = cgcolor('red')
  if photometrycolor eq 1 then photcolors = cgcolor(['CG11','CG2','CG3','CG4','CG5','CG6','CG7','CG8','CG9','CG10','CG1'])
  if r.name eq '2MASSJ12514562-8806328' then photcolors= cgcolor(['CG11','CG2','CG5'])
  syncolor = cgcolor('lime green')
  specolor = cgcolor('black')
  modelcolor = cgcolor('blue')
  symsize = 1.4

  !p.multi=[0,1,2]
  xmargin = !x.margin
  ymargin = !y.margin
  !x.margin = [9,2]
  !y.margin = [-4,1]

  downsampler = findgen(n_elements(lambda)) ;findgen(n_elements(lambda)/4)*4
  nm = name
  strreplace,nm,'TIC','TIC '

  if oc eq 1 then nameadd = '_OC' else nameadd = ''
  if long eq -1 then device,filename='Phot_fits/'+name+'_'+string(add,format="(I1)")+nameadd+'.eps',/encapsul,/color,xsize=30,ysize=10 else device,filename='Phot_fits/'+name+'_'+string(add,format="(I1)")+nameadd+'.eps',/encapsul,/color
  lamflam = cgsymbol('lambda') + 'F!L'+cgsymbol('lambda')+'!N'
  plot,lambda[downsampler],spec[downsampler],xrange=xrange,/xstyle,/ystyle,xtickn=strarr(10)+' ',charsize=charsize,charthick=charthick,$
  xthick=xthick,ythick=ythick,thick=thick,ytitle=lamflam+'(!Nerg cm!U-2!N s!U-1!N)',yrange=yrange,/ylog,/xlog,xticklen=0.05
  ;;xyouts,0.7,yrange[0]*2,nm,charsize=charsize-0.25,charthick=thick,alignment=0
  sharpcorners,thick=xthick,color=cgcolor('black')
  if add eq 1 then begin
     q = where(err gt median(err)*0.9995 and err lt median(err)*1.0005 and lambda gt 1.31 and lambda lt 1.45)
     if n_Elements(q) gt 10 then begin
        q = q[where((q-shift(q,1)) le 1)]
        oplot,lambda[q],(spec[q]),color=modelcolor,thick=thick
     endif
     q = where(err gt median(err)*0.9995 and err lt median(err)*1.0005 and lambda gt 1.75 and lambda lt 2.1)
     if n_Elements(q) gt 10 then begin
        q = q[where((q-shift(q,1)) le 1)]
        oplot,lambda[q],(spec[q]),color=modelcolor,thick=thick
     endif
     q = where(err gt median(err)*0.9995 and err lt median(err)*1.0005 and lambda lt 0.5)
     if n_Elements(q) gt 10 then begin
        q = q[where((q-shift(q,1)) le 1)]
        oplot,lambda[q],(spec[q]),color=modelcolor,thick=thick
     endif
  endif
  if add ge 1 then begin
     q = where(err gt median(err)*0.995 and err lt median(err)*1.005 and lambda gt 2.35)
     if n_Elements(q) gt 10 then oplot,lambda[q],(spec[q]),color=modelcolor,thick=thick
  endif

  oplot,lams[l],fluxes[l],color=syncolor,psym=8,symsize=symsize

  if photometrycolor eq 1 then begin
     syslist = (phot_system[sort(phot_system)])[uniq(phot_system[sort(phot_system)])]
     for jj = 0,n_elements(syslist)-1 do begin
        ll = where(phot_system eq syslist[jj])
        plotsym,0,thick=thick
        oploterror,lams[ll],phfluxes[ll],fwhms[ll],phfluxes_err[ll],color=photcolors[jj],psym=3,errthick=errthick,errcolor=photcolors[jj]
     endfor
  endif else oploterror,lams[l],phfluxes[l],fwhms[l],phfluxes_err[l],color=photcolor,psym=3,errthick=errthick,errcolor=photcolor
  plotsym,0,/fill



  if r.name eq 'EPIC204376071' then begin
     plotsym,1,/fill,thick=errthick
     oploterror,[22.375208],[5.0413475e-13],[2.0508831],[0],psym=8,color=photcolor,errthick=errthick,errcolor=photcolor,symsize=3.0
     plotsym,0,/fill
     oplot,[22.375208],[5.5032388e-14 ],psym=8,color=syncolor,symsize=symsize
  endif
  ;; values for EPIC204376071
  ;; template prediction = 5.5032388e-14
  ;; photometry = 5.0413475e-13

  if printname eq 1 then begin
     name = r.name
     strreplace,name,'_',' '
     xyouts,0.7,0.73,name,/normal,charsize=charsize,charthick=charthick
  endif
  if r.name eq 'HD_189733' then begin
     xyouts,3.2,0.5,'HD 189733',charsize=charsize,charthick=charthick
  endif
  if r.name eq 'HD_209458' then begin
     xyouts,3.0,0.5,'HD 209458',charsize=charsize,charthick=charthick
  endif
  xx_tmp = 13000
  yy_tmp = 9000
  if r.cns3 eq 'GJ 411' then begin
     xyouts,xx_tmp,yy_tmp,'Gl 411',charsize=charsize,charthick=charthick,/device
  endif
  if r.cns3 eq 'GJ 699' then begin
     xyouts,xx_tmp,yy_tmp,'Gl 699',charsize=charsize,charthick=charthick,/device
  endif
  if r.cns3 eq 'GJ 880' then begin
     xyouts,xx_tmp,yy_tmp,'Gl 880',charsize=charsize,charthick=charthick,/device
  endif
  if r.cns3 eq 'GJ 876' then begin
     xyouts,xx_tmp,yy_tmp,'Gl 876',charsize=charsize,charthick=charthick,/device
  endif

  sharpcorners,thick=xthick,color=cgcolor('black')

;;xyouts,2.5,median(spec[where(lambda gt 1 and lambda lt 3)]/(10.^factor)),string(rchisq_cover,format='(D5.2)'),charsize=charsize,alignment=0.5
  top = 1
  bottom = 0
  left = 0
  right = 1
  template = 'Template'
  ;;if stis eq 1 then op = 'NGSL' else op = 'SNIFS'
  ;;template = op+'+SpeX'
  ;;if soar eq 1 then template = 'SOAR+Template'
  if add ge 1 then begin
     ;;position = [2,1d-9]
     if r.name eq 'EPIC204376071' then position = [5,5d-13]
     if r.name eq 'HD222259B' then begin
        position = [2,5d-10]
        ttop = 0
        lleft = 1
        rright = 0
        bbottom = 1
     endif
     if r.name eq 'HD222259' then begin
        position = [2,1d-9]
        ttop = 0
        lleft = 1
        rright = 0
        bbottom = 1
     endif
     if r.name eq 'TYC7577-172-1' then begin
        position = [2,1d-9]
        ttop = 1
        lleft = 0
        rright = 1
        bbottom = 0
     endif
     if r.name eq 'PM_I20451-3120' then begin
        ttop = 1
        bbottom = 0
        lleft = 0
        rright = 1
     endif
     if r.name eq 'HD110082' then begin
        ttop = 1
        bbottom = 0
        lleft = 0
        rright = 1
     endif
     if r.name eq 'HD120411' then begin
        ttop = 1
        bbottom = 0
        lleft = 0
        rright = 1
     endif
     if r.name eq '2MASSJ16100501-2132318' then begin
        ttop = 0
        bbottom = 1
        lleft = 0
        rright = 1
     endif
     if r.name eq 'TYC3496-1082-1' then begin
        ttop = 0
        bbottom = 1
        lleft = 1
        rright = 0
     endif
     if photometrycolor eq 1 then begin
        if nolegend eq 0 then legend,[' ','              ','Photometry','Synthetic Photometry'],color=[specolor,modelcolor,photcolors[0],syncolor],psym=[0,0,1,8],charsize=charsize/1.25,box=0,textcolor=[specolor,modelcolor,photcolors[0],syncolor],top=0,right=0,left=1,bottom=1,thick=[0,0,thick,1]
        if nolegend eq 0 then legend,[template,'BT-SETTL Model',' ',' '],                       color=[specolor,modelcolor,photcolors[0],syncolor],psym=[0,0,1,8],charsize=charsize/1.25,box=0,textcolor=[specolor,modelcolor,photcolors[0],syncolor],top=0,right=0,left=1,bottom=1,thick=[10,10,0,0] ;;,top=top,right=right,left=left,bottom=bottom

     endif else begin
        ;;ttop = 0
        bbottom = 1
        lleft = 1
        ;;rright = 1
        if nolegend eq 0 then legend,[' ','              ','Photometry','Synthetic Photometry'],color=[specolor,modelcolor,photcolor,syncolor],psym=[0,0,8,8],charsize=charsize/1.25,box=0,textcolor=[specolor,modelcolor,photcolor,syncolor],top=ttop,right=rright,left=lleft,bottom=bbottom
        ;;ttop = 0
        bbottom = 1
        lleft = 1
        ;;rright = 1
        if nolegend eq 0 then legend,[template,'BT-SETTL Model',' ',' '],color=[specolor,modelcolor,photcolor,syncolor],psym=[0,0,8,8],charsize=charsize/1.25,box=0,thick=[10,10,0,0],textcolor=[specolor,modelcolor,photcolor,syncolor],top=ttop,right=rright,left=lleft,bottom=bbottom
        ;;ttop = 0
        bbottom = 1
        lleft = 1
        ;;rright = 1
     endelse
  endif else begin
    ;;stop
     if nolegend eq 0 then legend,['Spectrum','Photometry','Synthetic Photometry'],color=[specolor,photcolor,syncolor],psym=[0,8,8],thick=[10,0,0],top=ttop,right=rright,left=lleft,charsize=charsize/1.25,box=0,textcolor=[specolor,photcolor,syncolor]
  endelse

  if photometrycolor eq 1 then begin
     tmp = strtrim(syslist,2)
     ll = where(tmp eq 'SDSS') & if ll[0] ne -1 then tmp[ll] = 'SkyMapper'
     ll = where(tmp eq 'Tycho') & if ll[0] ne -1 then tmp[ll] = 'Tycho-2'
     ll = where(tmp eq 'Johnson') & if ll[0] ne -1 then tmp[ll] = 'APASS'
     ;tmp = tmp[where(tmp ne 'Wise')]
     if nolegend eq 0 then legend,tmp,psym=1+intarr(n_elements(tmp)),color=photcolors[0:n_elements(tmp)-1],charsize=charsize/1.25,box=0,textcolor=photcolors[0:n_elements(tmp)-1],thick=4+intarr(n_elements(tmp)),top=1,right=1,left=0,bottom=0

  endif

  if inset eq 1 then begin
     plot, lambda, spec/(10.^factor), xrange=[0.38,0.49], /ystyle,/xstyle, /noeras,ytickname=strarr(10)+' ',thick=thick,xthick=thick,ythick=thick,color=specolor,xtitle='Wavelength (!9m!3m)',charsize=0.9,position=[0.3,0.42,0.6,0.6] ;position=[0.62,0.52,0.94,0.77],
     sharpcorners,thick=xthick,color=specolor
  endif

  !y.margin=[3.5,4]
  sharpcorners,thick=xthick,color=cgcolor('black')
  xtickn = ['0.5',' ',' ',' ',' ','1','2',' ',' ','5',' ',' ',' ',' ','10']
  xtickv = [0.5,0.6,0.7,0.8,0.9,1.0,2.0,3,4,5,6,7,8,9,10]
  if pxrange[1] lt 10 then begin
     xtickn = ['0.5',' ',' ',' ',' ','1','2',' ',' ','5']
     xtickv = [0.5,0.6,0.7,0.8,0.9,1.0,2.0,3,4,5]
  endif
  if pxrange[1] lt 5 then begin
     xtickn = ['0.5',' ',' ',' ',' ','1','2']
     xtickv = [0.5,0.6,0.7,0.8,0.9,1.0,2.0]
  endif
  if r.name eq 'HD222259' then begin
     xtickn = ['0.5',' ',' ',' ',' ','1','2',' ',' ','5',' ',' ',' ',' ','10']
     xtickv = [0.5,0.6,0.7,0.8,0.9,1.0,2.0,3,4,5,6,7,8,9,10]
  endif
  if r.name eq 'EPIC204376071' then begin
     xtickn = ['0.5',' ',' ',' ',' ','1','2',' ',' ','5',' ',' ',' ',' ','10','20']
     xtickv = [0.5,0.6,0.7,0.8,0.9,1.0,2.0,3,4,5,6,7,8,9,10,20]
  endif
  if r.name eq 'PM_I20451-3120' then begin
     xtickn = ['0.5',' ','0.7',' ',' ','1','2',' ',' ','5']
     xtickv = [0.5,0.6,0.7,0.8,0.9,1.0,2.0,3,4,5]
  endif
  if oc eq 0 then begin                                   ;; sigmas
     yrange = [-3.,3.]
     if r.name eq 'HD222259' then yrange=[-3.5,3.5]
     if r.name eq 'HD222259B' then yrange=[-3.5,3.5]
     sigma = greek('sigma')+'!3' ;textoidl('\sigma')
     plot,[0],[0],xrange=xrange,yrange=yrange,charsize=charsize,charthick=charthick,xthick=5,ythick=5,xtitle='Wavelength (!9m!3m)',/xstyle,/ystyle,xticklen=0.06,ytitle='Residual ('+sigma+')',/xlog,xtickn=xtickn,xtickv=xtickv,xticks=n_elements(xtickv)-1,yminor=1
     if photometrycolor eq 1 then begin
        for jj = 0,n_elements(syslist)-1 do begin
           ll = where(phot_system eq syslist[jj])
           oplot,lams[ll],(phfluxes[ll]-(fluxes[ll]))/sqrt(fluxes_err[ll]^2.0+phfluxes_err[ll]^2.0),color=photcolors[jj],psym=errsym
        endfor
     endif else oplot,lams[l],(phfluxes[l]-(fluxes[l]))/sqrt(fluxes_err[l]^2.0+phfluxes_err[l]^2.0),psym=errsym

  endif else begin
     magdiff = -2.5*alog10(fluxes/phfluxes)
     nmonte = 1d5
     magdiff_err = dblarr(n_elements(magdiff))
     for jjj = 0,n_elements(phfluxes)-1 do begin
        tmp1 = phfluxes[jjj]+phfluxes_err[jjj]*randomn(seed,nmonte)
        tmp2 = fluxes[jjj]+fluxes_err[jjj]*randomn(seed,nmonte)
        tmp3 = -2.5*alog10(tmp1/tmp2)
        magdiff_err[jjj] = robust_sigma(tmp3)
     endfor
     q = where(lams gt xrange[0] and lams lt xrange[1])
     min = min([magdiff[q]-magdiff_err[q]])
     max = max([magdiff[q]+magdiff_Err[q]])
     yrange = [-1d0*max(abs([min,max])),max(abs([min,max]))]
     ploterror,lams[l],magdiff[l],lams[l]*0d0,magdiff_err[l],psym=errsym,errthick=thick,xrange=xrange,/xstyle,/ystyle,xtitle='Wavelength (!9m!3m)',charsize=charsize,charthick=charthick,ytitle='Residual (mag)',/xlog,xtickn=xtickn,xtickv=xtickv,xticks=n_elements(xtickv)-1,yrange=yrange,xthick=5,ythick=5,xticklen=0.05

  endelse
  oplot,[0.01,100],[0,0],linestyle=2,thick=5,color=cgcolor('black')
  oplot,[0.01,100],[0,0],linestyle=2,thick=5,color=cgcolor('black')
  sharpcorners,thick=xthick,color=cgcolor('black')
  device,/close

  !x.margin = xmargin
  !y.margin = ymargin
  set_plot,'x'

  fac = 10000d0
  spec/=lambda*fac
  ;;sp_err/=lambda*fac
  phfluxes/=lams*fac
  phfluxes_err/=lams*fac
  FLUXES/=lams*fac
  FLUXES_err/=lams*fac

END
