pro gstruc_add,tstr
;;  This adds new elements to gstruct

ntstr = n_elements(tstr)

gstruc_schema = {lon:999999.,lat:999999.,rms:999999.,noise:999999.,$
                 par:fltarr(3)+999999.,sigpar:fltarr(3)+999999.,glon:999999.,glat:999999.}

DEFSYSV,'!gstruc',exists=gstruc_exists
if not keyword_set(gstruc_exists) then $
  DEFSYSV,'!gstruc',{data:replicate(gstruc_schema,100000L),ndata:100000LL,count:0LL,revindex:lon64arr(100000L)-1,$
                     lonstart:fltarr(100000L)+999999.,latstart:fltarr(100000L)+999999.,indstart:lon64arr(100000L)-1,ngauss:lonarr(100000L)-1,pcount:0LL}
;; data: large data structure
;; ndata: number of elements of data
;  count:  the number of elements in DATA were using currently, also
;           the index of the next one to stat with.
;; revindex: reverse index, takes you from data index to lonstart/latstart/ngauss
;; lonstart/latstart/ngauss: give the lon/lat value at the start of a
;;            sequence of Gaussians, and the number of gaussians
;; indstart: the index in DATA where the gaussians for this position start
;; pcount: the number of spatial positions

;; Add new elements
if ntstr+!gstruc.count gt !gstruc.ndata then begin
  print,'Adding more elements to GSTRUC'
  ;; Old structures/arrays
  data = !gstruc.data
  ndata = !gstruc.ndata
  count = !gstruc.count
  revindex = !gstruc.revindex
  lonstart = !gstruc.lonstart
  latstart = !gstruc.latstart
  indstart = !gstruc.indstart
  pcount = !gstruc.pcount
  ;; Make new structures/arrays
  new_data = replicate(gstruc_schema,ndata+100000L)
  new_ndata = long64(n_elements(new_data))
  new_revindex = lon64arr(ndata+100000L)-1
  new_lonstart = fltarr(ndata+100000L)+999999.
  new_latstart = fltarr(ndata+100000L)+999999.
  new_indstart = lon64arr(ndata+100000L)-1
  new_ngauss = lonarr(ndata+100000L)-1
  ;; Stuff in the old values
  new_data[0:ndata-1] = data
  new_revindex[0:ndata-1] = revindex
  new_lonstart[0:pcount-1] = lonstart
  new_latstart[0:pcount-1] = latstart
  new_indstart[0:pcount-1] = indstart
  new_ngauss[0:pcount-1] = ngauss
  DEFSYSV,'!gstruc',{data:new_data,ndata:new_ndata,count:count,revindex:new_revindex,lonstart:new_lonstart,$
                     latstart:new_latstart,indstart:new_indstart,ngauss:new_ngauss,pcount:pcount}
endif
;; Stuff in the new data
!gstruc.data[!gstruc.count:!gstruc.count+ntstr-1] = tstr
!gstruc.revindex[!gstruc.count:!gstruc.count+ntstr-1] = !gstruc.pcount
!gstruc.lonstart[!gstruc.pcount] = tstr[0].lon
!gstruc.lonstart[!gstruc.pcount] = tstr[0].lat
!gstruc.indstart[!gstruc.pcount] = !gstruc.count
!gstruc.ngauss[!gstruc.pcount] = ntstr
!gstruc.count += ntstr
!gstruc.pcount += 1

end
