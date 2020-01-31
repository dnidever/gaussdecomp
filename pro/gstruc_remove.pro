pro gstruc_remove,ind

;; This "removes" indices from gstruc by shifting everything over
;;  this assumes that the IND indices is contiguous

nind = n_elements(ind)
rem0 = min(ind)
rem1 = max(ind)

;; UPDATE DATA
;; At the end, nothing to shift, just zero-out
if rem1 eq !gstruc.count-1 then begin
  ;; Zero out the removed elements
  ntags = n_tags(!gstruc.data)
  for i=0,ntags-1 do !gstruc.data[rem0:rem1].(i)=999999.

;; Shift everything above REM1 to start at REM0
endif else begin
  above0 = rem1+1
  above1 = !gstruc.count-1
  nabove = above1-above0+1
  temp = !gstruc.data[above0:above1]
  ;; Zero out the old positions
  ntags = n_tags(gstruc)
  for i=0,ntags-1 do !gstruc.data[above0:above1].(i)=999999.
  ;; Stuff the data in the new positions
  !gstruc.data[rem0:rem0+nabove-1] = temp
endelse

;; UPDATE LONSTART/LATSTART/NGAUSS arrays
;; lonstart, latstart, indstart, ngauss, pcount
;; Use REVINDEX
pind = !gstruc.revindex[ind[0]]
;; At the end, nothing to shift, just zero-out
if pind eq !gstruc.pcount-1 then begin
  ;; Zero out the remaining elements
  !gstruc.lonstart[pind] = 999999.
  !gstruc.latstart[pind] = 999999.
  !gstruc.indstart[pind] = -1
  !gstruc.ngauss[pind] = -1

;; Shift everything above PREM1 to start at PREM0
endif else begin
  pabove0 = pind+1
  pabove1 = !gstruc.pcount-1
  npabove = pabove1-pabove0+1
  ;; LONSTART
  temp_lonstart = !gstruc.lonstart[pabove0:pabove1]
  !gstruc.lonstart[pabove0:pabove1] = 999999.
  !gstruc.lonstart[pind:pind+npabove-1] = temp_lonstart
  ;; LATSTART
  temp_latstart = !gstruc.latstart[pabove0:pabove1]
  !gstruc.latstart[pabove0:pabove1] = 999999.
  !gstruc.latstart[pind:pind+npabove-1] = temp_latstart
  ;; INDSTART
  temp_indstart = !gstruc.indstart[pabove0:pabove1]
  !gstruc.indstart[pabove0:pabove1] = -1
  !gstruc.indstart[pind:pind+npabove-1] = temp_indstart
   ;; NGAUSS
  temp_ngauss = !gstruc.ngauss[pabove0:pabove1]
  !gstruc.ngauss[pabove0:pabove1] = -1
  !gstruc.ngauss[pind:pind+npabove-1] = temp_ngauss
  ;; Renumber NGAUSS
  ;;  convert the old index to the new index
  ;;  newindex = converter[oldindex]
  indtokeep = l64indgen(!gstruc.count)
  REMOVE,ind,indtokeep
  index_converter = lon64arr(!gstruc.count)-1
  index_converter[indtokeep] = l64indgen(!gstruc.count-nind)
  ;; only convert the ones that were moved
  moved_old_ngauss_values = !gstruc.ngauss[pind:pind+npabove-1]
  !gstruc.ngauss[pind:pind+npabove-1] = index_convert[moved_old_ngauss_values]
endelse

;; UPDATE REVINDEX
;;  shift them first

;; At the end, nothing to shift, just zero-out
if rem1 eq !gstruc.count-1 then begin
  !gstruc.revindex[ind] = -1  ; zero out the removed elements

;; Shift everything above REM1 to start at REM0
endif else begin
  above0 = rem1+1
  above1 = !gstruc.count-1
  nabove = above1-above0+1
  temp = !gstruc.revindex[above0:above1]
  ;; Zero out the old positions
  !gstruc.revindex[above0:above1] = -1
  ;; Stuff the data in the new positions
  !gstruc.revindex[rem0:rem0+nabove-1] = temp
  ;; Renumber REVINDEX
  ;;  convert the old index to the new index
  ;;  newindex = convert[oldindex]
  indtokeep = l64indgen(!gstruc.pcount)
  REMOVE,pind,indtokeep
  index_converter = lon64arr(!gstruc.pcount)-1
  index_converter[indtokeep] = l64indgen(!gstruc.pcount-1)
  ;; only convert the ones that were moved
  moved_old_revindex_values = !gstruc.revindex[rem0:rem0+nabove-1]
  !gstruc.revindex[rem0:rem0+nabove-1] = index_convert[moved_old_revindex_values]
endelse

;; Reduce counters
!gstruc.pcount--
!gstruc.count -= nind

end
