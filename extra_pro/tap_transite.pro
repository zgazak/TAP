pro tap_transite, times, x, flux, trend=trend			  
  ;; Computes a transit lightcurve, normalized to unity, as a
  ;; function of time, t, usually in HJD or BJD.
  ;;
  ;; Input parameters (x) are:
  ;; x(0) = per = period in days
  ;; x(1) = inc = inclination in degrees
  ;; x(2) = rsource = R_*/a 
  ;; x(3) = p = R_p/R_*
  ;; x(4) = t0 = mid-point of transit (JD)
  ;; x(5) = u1
  ;; x(6) = u2
  ;; x(7) = ecc
  ;; x(8) = pomega [RADIANS!]
  
  ;; z0=sqrt(x(1)^2+((t-x(3))*x(0))^2)	
  per = x[0]
  inc = x[1]*(!dpi/180d0)
  
  rsource = x[2]
  p = x[3]
;;; t0 (x[4]) comes in as T_mid, so I need to convert to Tperi for this code
  u1 = x[5]
  u2 = x[6]
  ecc = x[7]
  pom = x[8]

  if ecc gt 0 then begin
     f = !dpi/2d0-pom
     E = 2d0*atan(sqrt((1d0-ecc)/(1d0+ecc))*sin(f/2d0), cos(f/2d0))
     n = 2d0*!dpi/per
     tperi = x[4] - (E-ecc*sin(E))/n
  endif else tperi = x[4] + per/4d0
  
  m=(2.d0*!dpi*(times-tperi)/per) mod (2.d0*!dpi) ; (2.39)
  
  if(ecc ne 0d0) then begin
     kepler_ma,m,ecc,f          ; solved (2.64)
  endif else begin
     f=m
  endelse
  radius=(1d0-ecc^2)/(1d0+ecc*cos(f))           ; (2.20)
  gmsun = 1.32712440018d26                    ; cm^3/s^2
  rpsky=radius*sqrt(1d0-(sin(inc)*sin(pom+f))^2) ; from (2.122) 
  
  occultquad_vec,rpsky/rsource,u1,u2,p,flux,F0
  
  if keyword_set(trend) then begin
     th = (times-min(times))*24d0
     flux /= poly(th, [x[9],x[10]])
   ;  F0 += poly(th, [x[9]-1d0,x[10]])
  endif
end
