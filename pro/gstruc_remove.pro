pro gstruc_remove,ind

;; This "removes" indices from gstruc by shifting everything over
;;  this assumes that the IND indices is contiguous

nind = n_elements(ind)
rem0 = min(ind)
rem1 = max(ind)

;; At the end, nothing to shift, just zero-out
if rem1 eq !gstruc.count-1 then begin
  ;; Zero out the removed elements
  ntags = n_tags(!gstruc.data)
  for i=0,ntags-1 do !gstruc.data[rem0:rem1].(i)=999999.

;; Shift everything above REM1 to star at REM0
endif else begin
  above0 = rem1+1
  above1 = gcount-1
  nabove = above1-above0+1
  !gstruc.data[rem0:rem0+nabove-1] = !gstruc.data[above0:above1]
  ;; Zero out the extra elements
  ntags = n_tags(gstruc)
  for i=0,ntags-1 do !gstruc.data[above0:above1].(i)=999999.
endelse

;; shift the lonstart/latstart/ngauss arrays
;; lonstart, latstart, indstart, ngauss, pcount
begind = where(!gstruc.indstart eq min(ind),nbegind)
if nbegind eq 0 then stop,'Cannot find the beginning index for this position'


!gstruc.count -= nind

end
