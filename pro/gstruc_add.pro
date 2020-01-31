pro gstruc_add,tstr
;;  This adds new elements to gstruct

ntstr = n_elements(tstr)

gstruc_schema = {lon:999999.,lat:999999.,rms:999999.,noise:999999.,$
                 par:fltarr(3)+999999.,sigpar:fltarr(3)+999999.,glon:999999.,glat:999999.}

DEFSYSV,'!gstruc',exists=gstruc_exists
if not keyword_set(gstruc_exists) then $
  DEFSYSV,'!gstruc',{data:replicate(gstruc_schema,100000L),ndata:100000LL,count:0LL}

;; Add new elements
if ntstr+!gstruc.count gt !gstruc.ndata then begin
  print,'Adding more elements to GSTRUC'
  data = !gstruc.data
  ndata = !gstruc.ndata
  cont = !gstruc.count
  newdata = replicate(gstruc_schema,ndata+100000L)
  newdata[0:ndata-1] = data
  nnewdata = long64(n_elements(newdata))
  DEFSYSV,'!gstruc',{data:newdata,ndata:nnewdata,count:count}
endif
;; Stuff in the new data
!gstruc.data[!gstruc.count:!gstruc.count+ntstr-1] = tstr
!gstruc.count += ntstr

end
