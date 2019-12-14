pro hinoise,v,spec,noise

;  This program estimates the noise in an HI spectrum
;
;  INPUT
;   v      Array of velocity
;   spec   Array of spectrum
;
;  OUTPUT
;   noise  The noise level of the spectrum
;
; When there is a problem in this program it returns:
;  noise = 999999.
; 
; Created by David Nidever April 2005

nv = n_elements(v)
nspec = n_elements(spec)

; Bad Input Values
if (n_params() eq 0) or (nv ne nspec) then begin
  print,'Syntax - hinoise,v,spec,noise'
  noise = 999999.
  return
endif

vmax = max(v)
vmin = min(v)

; Velocity array not long enough
;if (vmin gt -250.) or (vmax lt 250.) then begin
;if (vmin gt -30.) or (vmax lt 30.) then begin
;if (vmin gt -250.) or (vmax lt 100.) then begin
if (vmin gt -100.) or (vmax lt 100.) then begin
  print,'The velocity range must go beyond [-100,100] km/s'
  noise = 999999.
  return
endif

; Smoothing the data
smspec2 = savgolsm(spec, [2,2,2])
smspec4 = savgolsm(spec, [4,4,2])
smspec16 = savgolsm(spec, [16,16,2])

; Points at high velocity
;gd1 = where(abs(v) gt 250.,ngd1)
;gd1 = where(abs(v) gt 30.,ngd1)
;gd1 = where(v gt 100 or v lt -250.,ngd1)
gd1 = where(v gt 200 or v lt -200.,ngd1)
diff = (spec-smspec16)(gd1)

; Points below the threshold
thresh = 5.*stdev(diff) > 0.5
gd2 = where(diff lt thresh, ngd2)
;if ngd2 eq 0 then stop

; Estimating the noise
;noise = stddev(diff(gd2))
noise = MAD(diff[gd2])

;stop

end
