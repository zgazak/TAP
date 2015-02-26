pro solveredwv,d,s_r,s_w,redcomp=redcomp,whitecomp=whitecomp,alpha=alpha,sol=sol,silent=silent,zeropad=zeropad

;+
; NAME:
;       SOLVEREDWV
; PURPOSE:
;       Given a data vector, D, this algorithm attempts to find the best fit model
;       describing the data as 1/f noise plus white noise
; EXPLANATION:
;       Uses the wavelet technique as described by Carter & Winn
;       (2009) to determine parameters of noise formed as an additive
;       combination of (Gaussian) noise with power spectral density proportional to 1/f and
;       (Gaussian) white noise.  In particular, AMOEBA is used to maximize a
;       likelihood that is a function of two parameters sigma_r and
;       sigma_w which are related to the standard deviations of the 1/f
;       and white components, respectively.
;
;       Additionally, a wavelet filter is applied to separate the two
;       components.  These results may be returned to the caller as optional outputs.
;
; CALLING SEQUENCE:
;       solveredwv,d, s_r, s_w,[redcomp=redcomp, whitecomp=whitecomp,
;       alpha=alpha, sol=sol, /silent]
;
; INPUTS:
;       D       - Data vector, length must be a power of two
;       SIGMA_R - Red noise amplitude (not equal to 1/f-component
;                 RMS): If 2-element vector, then this input gives the
;                 approximate range for SIGMA_R in which the solution
;                 lies.  If scalar, SIGMA_R is fixed to its input
;                 value (i.e., it does not vary in AMOEBA).
;       SIGMA_W - White noise amplitude parameter (approximately equal
;                 to white-component RMS): If 2-element vector, then this input
;                 input gives the approximate range for SIGMA_W in 
;                 which the solution lies.  If scalar, SIGMA_W is
;                 fixed to its input value (i.e., it does not vary in
;                 AMOEBA).
; OPTIONAL INPUTS:
;       SILENT - Suppresses output.
;
; OUTPUTS:
;       Prints the best fit values for SIGMA_R and SIGMA_W.  Also
;       reports the RMS ratio of the 1/f component to the white
;       component [referred to as "alpha" by Carter & Winn (2009)]
;
; OPTIONAL OUTPUTS:
;       REDCOMP - The best fit 1/f-component solution for data vector d (of the
;                 same length as d).
;       WHITECOMP - The best fit white-component solution for data
;                 vector d (of the same length as d).
;       SOL     - 2-element vector of the best fit SIGMA_R (SOL[0]) and
;                 SIGMA_W (SOL[1])
;       ALPHA   - Number giving the RMS ratio of the 1/f component to
;                 the white component.
;
; EXAMPLES:
;
;       EXAMPLE 1: White noise.
;
;       IDL> x = randomn(seed,1024)
;       IDL> solveredwv,x,[0.0,1.0],[0.5,1.5],sol=sol,alpha=alpha,$
;            redcomp=redcomp,whitecomp=whitecomp
;       Sigma_r = -0.032677
;       Sigma_w = 0.994064
;       Gamma = 1.000000
;       RMS of white component is 0.993101
;       RMS of red component is 0.000021
;       Ratio of red to white RMS is 0.000021
;       IDL> 
;
;       EXAMPLE 2: Correlated noise (RMS = 0.40) + white noise (RMS = 1.0)
;
;       IDL> X = smooth(randomn(seed,1024),4)+randomn(seed,1024)
;       IDL> solveredwv,x,[0.0,1.0],[0.5,1.5],sol=sol,alpha=alpha,$
;            redcomp=redcomp,whitecomp=whitecomp
;       Sigma_r = 6.238085
;       Sigma_w = 1.000332
;       Gamma = 1.000000
;       RMS of white component is 0.951405
;       RMS of red component is 0.215253
;       Ratio of red to white RMS is 0.226248
;
; NOTES:
;       Refer to the paper by Carter & Winn (2009) for theory and
;       details.  In the notation of that paper, gamma=1 for this
;       algorithm.
;
; REVISION HISTORY:
;       Written,    J. Carter               September, 2009
;       Made default zero-padding option    November, 2009
;       Corrected non-short-circuited ifs!      
;               Thanks to E. Ford           January, 2010
;-   

  common solvewv, x,var,sigma_r,sigma_w,gamma

  sigma_r = s_r
  sigma_w = s_w
  gamma = 1

  els = n_elements(d)
  pow = ceil(alog(els)/alog(2.))
  if (2^pow ne els and keyword_set(zeropad)) then begin
     diff = 2^pow-els
     left = floor(diff/2.)
     right = diff-left
     x = [replicate(0d,left),d,replicate(0d,right)]
  endif else x = d

  var = replicate(0,3)

  if (n_elements(sigma_r) ne 1) then var(0) = 1
  if (n_elements(sigma_w) ne 1) then var(1) = 1
  if (n_elements(gamma) ne 1) then var(2) = 1

  k = 0
  par = dindgen(total(var))
  scale = dindgen(total(var))
  if (var(0) eq 1) then begin
     par(k) = mean(sigma_r)
     scale(k++) = sigma_r(1)-sigma_r(0)
  endif
  if (var(1) eq 1) then begin
     par(k) = mean(sigma_w)
     scale(k++) = sigma_w(1)-sigma_w(0)
  endif
  if (var(2) eq 1) then begin
     par(k) = mean(gamma)
     scale(k++) = gamma(1)-gamma(0)
  endif

  
  result = amoeba(1.0e-5,function_name='minfuncwv',p0=par,scale=scale)
  redcomp = filterredwv(x,sigma_r,sigma_w)
  if keyword_set(zeropad) then begin
     if left ne 0 and right ne 0 then begin
        x = d
        ;whitecomp = whitecomp(left:(els+left))
        redcomp = redcomp(left:(els+left))
     endif
  endif
  whitecomp = x - redcomp

  if (not(keyword_set(silent))) then begin
     print,sigma_r,format="('Sigma_r = ',F0)"
     print,sigma_w,format="('Sigma_w = ',F0)"
     print,gamma,format="('Gamma = ',F0)"

     print,stddev(whitecomp),format="('RMS of white component is ', F0)"
     print,stddev(redcomp),format="('RMS of red component is ',F0)"

     print,stddev(redcomp)/stddev(whitecomp),format="('Ratio of red to white RMS is ',F0)"
  endif

  sol = [sigma_r,sigma_w]
  alpha = stddev(redcomp)/stddev(whitecomp)

end


  
  
