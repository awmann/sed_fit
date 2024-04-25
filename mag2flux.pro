
function mag2flux, band,mag,err
;sdss softening parameters
;;b_u = 1.4D-10
  band = strtrim(band,2)
  case band of
     'u': b = 1.4D-10
     'g': b = 0.9D-10
     'r': b = 1.2D-10
     'i': b = 1.8D-10
     'z': b = 7.4D-10
  endcase

  jansk = 3631.*2.*b*sinh((mag/(-2.5/alog(10.))) - alog(b))
  nmonte = 100
  newmags = mag+err*randomn(seed,nmonte)
  jansks = 3631.*2.*b*sinh((newmags/(-2.5/alog(10.))) - alog(b))

;calculate asinh AB mags
  ;;gmag = -1.*(2.5/alog(10))*[asinh((janskys[0]/3631.)/(2*b_g))+alog(b_g)]
  ;;rmag = -1.*(2.5/alog(10))*[asinh((janskys[1]/3631.)/(2*b_r))+alog(b_r)]
  ;;imag = -1.*(2.5/alog(10))*[asinh((janskys[2]/3631.)/(2*b_i))+alog(b_i)]
  ;;zmag = -1.*(2.5/alog(10))*[asinh((janskys[3]/3631.)/(2*b_z))+alog(b_z)]
  ;;mags=[gmag,rmag,imag,zmag]g

  return, [jansk,stdev(jansks)]
end
