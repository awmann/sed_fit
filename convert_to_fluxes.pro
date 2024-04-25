PRO convert_to_fluxes,struct,system,band,mag,mag_err,flux,flux_err,leff
  COMPILE_OPT idl2, HIDDEN
  if system eq 'Eggen' and band eq 'V' then system = 'Johnson'
  if system eq 'Cousins' and band eq 'U' then system = 'Johnson'
  if system eq 'Cousins' and band eq 'B' then system = 'Johnson'
  if system eq 'Cousins' and band eq 'V' then system = 'Johnson'
  if system eq 'Straizys' then system = 'Vilnius'
  if system eq 'Cape' and (band eq 'V' or band eq 'B') then system = 'Johnson'
  if system eq 'Johnson' and band eq 'I' then begin ;; this system doesnt work
     mag = -99
     mag_err = -99
  endif else begin
     if system eq 'SDSS' then begin
        mag_err*=2.0
        mag_err = sqrt(mag_err^2.0+0.0175^2.0)
        l = where(STRLOWCASE(struct.system) eq STRLOWCASE(system) and STRLOWCASE(struct.band) eq STRLOWCASE(band))
        tmp = mag2flux(band,mag,mag_err)
        flux = tmp[0]*2.99792458d-5/(struct[l].leff*1d4)^2.
        flux_err = tmp[1]*2.99792458d-5/(struct[l].leff*1d4)^2.
        leff = struct[l].leff
     endif else begin
        if STRLOWCASE(strtrim(system,2)) eq '2mass' then begin
           if band eq 'J' then adder = 0.001
           if band eq 'H' then adder = -0.019
           if band eq 'K' then adder = 0.017
        endif else adder = 0.
        if system eq 'Hipparcos' then begin
           mag_err = sqrt(mag_err^2.0+0.01^2.0)
        endif
        if system eq 'Tycho' then begin
           mag_err = sqrt(mag_err^2.0+0.01^2.0)
        endif
        if system eq 'Wise' then begin
           if band eq 'W1' then mag_err = sqrt(mag_err^2.0+0.06^2.0)
           if band eq 'W2' then mag_err = sqrt(mag_err^2.0+0.10^2.0)
           if band eq 'W3' then mag_err = sqrt(mag_err^2.0+0.04^2.0)
           if band eq 'W4' then mag_err = sqrt(mag_err^2.0+0.07^2.0)
        endif        
        if band eq 'L' then begin
           mag_err = 0.1
        endif
        l = where(STRLOWCASE(struct.system) eq STRLOWCASE(system) and STRLOWCASE(struct.band) eq STRLOWCASE(band))
        if l[0] eq -1 then begin
           flux = -99
           flux_err = -99
        endif else begin
           if system eq 'Gaia' then mag_err+=0.025 ;; ZP errors
           flux = struct[l].zp*10d0^(-0.4*(mag+adder))
           flux_pl = struct[l].zp*10d0^(-0.4*(mag+adder+mag_err))
           flux_mi = struct[l].zp*10d0^(-0.4*(mag+adder-mag_err))
           flux_err = abs(flux_pl-flux_mi)/2.0
           leff = struct[l].leff
           if finite(flux_err) eq 0 then stop
           if system eq 'Wise' then begin ;; corrections
              if band eq 'W1' then begin
                 ;;flux/=1.0315
                 ;;flux_pl/=1.0315
                 ;;flux_mi/=1.0315
              endif
              if band eq 'W2' then begin
                 ;;flux/=1.1041
                 ;;flux_pl/=1.1041
                 ;;flux_mi/=1.1041
              endif
              if band eq 'W3' then begin
                 flux/=0.8906
                 flux_pl/=0.8906
                 flux_mi/=0.8906
              endif
              if band eq 'W4' then begin
                 flux/=0.7935
                 flux_pl/=0.7935
                 flux_mi/=0.7935
              endif
           endif
           if system eq 'Eggen' and (band eq 'R' or band eq 'I') then flux/=0.9667
        endelse
     endelse
  endelse
;;if leff gt 0.54 and leff lt 0.55 then stop
END
