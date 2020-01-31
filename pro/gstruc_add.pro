pro gstruc_add,tstr
;;  This adds new elements to gstruct

ntstr = n_elements(tstr)

gstruc_schema = {lon:999999.,lat:999999.,rms:999999.,noise:999999.,$
                 par:fltarr(3)+999999.,sigpar:fltarr(3)+999999.,glon:999999.,glat:999999.}

DEFSYSV,'!gstruc',exists=gstruc_exists
if not keyword_set(gstruc_exists) then $
  DEFSYSV,'!gstruc',{data:replicate(gstruc_schema,100000L),ndata:100000LL,count:0LL,$
                     lonstart:fltarr(100000L),latstart:fltarr(100000L),indstart:lon64arr(100000L),ngauss:lonarr(100000L),pcount:0LL}
;; data: large data structure
;; ndata: number of elements of data
;  count:  the number of elements in DATA were using currently, also
;           the index of the next one to stat with.
;; lonstart/latstart/ngauss: give the lon/lat value at the start of a
;;            sequence of Gaussians, and the number of gaussians
;; indstart: the index in DATA where the gaussians for this position start
;; pcount: the number of spatial positions

;; Add new elements
if ntstr+!gstruc.count gt !gstruc.ndata then begin
  print,'Adding more elements to GSTRUC'
  data = !gstruc.data
  ndata = !gstruc.ndata
  cont = !gstruc.count
  lonstart = !gstruc.lonstart
  latstart = !gstruc.latstart
  indstart = !gstruc.indstart
  pcount = !gstruc.pcount
  newdata = replicate(gstruc_schema,ndata+100000L)
  newdata[0:ndata-1] = data
  nnewdata = long64(n_elements(newdata))
  newlonstart = fltarr(ndata+100000L)
  newlatstart = fltarr(ndata+100000L)
  newindstart = lon64arr(ndata+100000L)
  newngauss = lonarr(ndata+100000L)
  newlonstart[0:pcount-1] = lonstart
  newlatstart[0:pcount-1] = latstart
  newindstart[0:pcount-1] = indstart
  newngauss[0:pcount-1] = ngauss
  DEFSYSV,'!gstruc',{data:newdata,ndata:nnewdata,count:count,lonstart:lonstart,latstart:latstart,indstart:indstart,ngauss:ngauss,pcount:pcount}
endif
;; Stuff in the new data
!gstruc.data[!gstruc.count:!gstruc.count+ntstr-1] = tstr
!gstruc.lonstart[!gstruc.pcount] = tstr[0].lon
!gstruc.lonstart[!gstruc.pcount] = tstr[0].lat
!gstruc.indstart[!gstruc.pcount] = !gstruc.count
!gstruc.ngauss[!gstruc.pcount] = ntstr
!gstruc.count += ntstr
!gstruc.pcount += 1

end
