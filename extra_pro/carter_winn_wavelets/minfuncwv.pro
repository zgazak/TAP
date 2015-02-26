function minfuncwv,p


;+
; NAME:
;       MINFUNCWV
; PURPOSE:
;       Helper function for SOLVEREDWV, not meant to be called separately.
;
; NOTES:
;       Refer to the paper by Carter & Winn (2009) for theory and
;       details.  In the notation of that paper, gamma=1 for this
;       algorithm.
;
; REVISION HISTORY:
;       Written,    J. Carter               September, 2009
;-   

   common solvewv, x,var,sigma_r,sigma_w,gamma

   i = 0
   sigma_r = (var(0) eq 1 ? p(i++) : sigma_r)
   sigma_w = (var(1) eq 1 ? p(i++) : sigma_w)
   gamma   = (var(2) eq 1 ? p(i++) : gamma)

   return, -1*waveletlike(x,sigma_r,sigma_w)

end
