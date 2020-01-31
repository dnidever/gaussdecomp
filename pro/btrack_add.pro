pro btrack_add,track
;;  This adds new elements to btrack

ntrack = n_elements(track)

np = 99
btrack_schema = {count:999999.,lon:999999.,lat:999999.,rms:999999.,noise:999999.,par:fltarr(np)+999999,$
                 guesspar:fltarr(np)+999999.,guesslon:999999.,guesslat:999999.,back:999999.,redo:999999.,$
                 redo_fail:999999.,skip:999999.,lastlon:999999.,lastlat:999999.}

DEFSYSV,'!btrack',exists=btrack_exists
if not keyword_set(btrack_exists) then $
  DEFSYSV,'!btrack',{data:replicate(btrack_schema,100000L),ndata:100000LL,count:0LL}

;; Add more gaussian parameter elements, PAR
if n_elements(track[0].par) gt n_elements(!btrack.data[0].par) then begin
  print,'Increasing the number elements of PAR in BTRACK'
  data = !btrack.data
  ndata = !btrack.ndata
  count = !btrack.count
  gbtrack,data,track  ;; this fixes the structure
  ;; New structure
  np = n_elements(track[0].par)
  btrack_schema = {count:999999.,lon:999999.,lat:999999.,rms:999999.,noise:999999.,par:fltarr(np)+999999,$
                   guesspar:fltarr(np)+999999.,guesslon:999999.,guesslat:999999.,back:999999.,redo:999999.,$
                   redo_fail:999999.,skip:999999.,lastlon:999999.,lastlat:999999.}  
  DEFSYSV,'!btrack',{data:data,ndata:ndata,count:count}
endif

;; Add new elements
if ntrack+!btrack.count gt !btrack.ndata then begin
  print,'Adding more rows to BTRACK'
  data = !btrack.data
  ndata = !btrack.ndata
  cont = !btrack.count
  newdata = replicate(btrack_schema,ndata+100000L)
  newdata[0:ndata-1] = data
  nnewdata = long64(n_elements(newdata))
  DEFSYSV,'!btrack',{data:newdata,ndata:nnewdata,count:count}
endif
;; Stuff in the new data
!btrack.data[!btrack.count:!btrack.count+ntrack-1] = track
!btrack.count += ntrack

end
