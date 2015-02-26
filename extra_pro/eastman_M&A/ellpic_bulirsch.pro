function ellpic_bulirsch,n,k
; Computes the complete elliptical integral of the third kind using
; the algorithm of Bulirsch (1965):
kc=sqrt(1d0-k^2) &  p=n+1d0
if(min(p) lt 0d0) then print,'Negative p'
m0=1d0 & c=1d0 & p=sqrt(p) & d=1d0/p & e=kc

while 1 do begin
    f = c & c = d/p+c & g = e/p & d = 2d0*(f*g+d)
    p = g + p & g = m0 & m0 = kc + m0
    if(max(abs(1d0-kc/g)) gt 1.d-8) then begin
        kc = 2d0*sqrt(e) & e=kc*m0
    endif else return,0.5d0*!dpi*(c*m0+d)/(m0*(m0+p))
endwhile

end
