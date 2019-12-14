pro gass_gauss_comb

; This combines the Gaussians for the 4 GASS subcubes.


;gass1.in on stream
;-------------------
restore,'gass1.dat'
add_tag,gstruc,'subcube',0L,gstruc
gstruc.subcube = 1
add_tag,gstruc,'X',0L,gstruc
add_tag,gstruc,'Y',0L,gstruc
gstruc.y = gstruc.lon+1
gstruc.x = gstruc.lat+2
newgstruc = gstruc
add_tag,gstruc,'MLON',0.0
add_tag,gstruc,'MLAT',0.0
gstruc.mlon = gstruc.glon
gstruc.mlat = gstruc.glat
mag2gal,gstruc.mlon,gstruc.mlat,glon,glat
gstruc.glon = glon
gstruc.glat = glat
; This cube spans -92.0 < MLON < -48.0
; Use -92.0 to -50.0
gd = where(gstruc.mlon le -50.0,ngd)
newgstruc = gstruc[gd]

;; Now write the ASCII file
;id = lindgen(n_elements(gstruc))+1
;fmt = '(I6,I5,I5,2F9.3,2F8.3,F9.3,F9.3,F9.3)' 
;;fmt = '(I8,2F7.1,2F8.3,F9.3,F9.3,F9.3)' 
;file = 'gass1.txt'
;writecol,file,id,gstruc.x,gstruc.y,gstruc.glon,gstruc.glat,gstruc.rms,gstruc.noise,$
;  gstruc.par[0],gstruc.par[1],gstruc.par[2],fmt=fmt
;
;undefine,header
;push,header,'#   NUM   X    Y   GLON      GLAT     RMS    NOISE   HEIGHT    VLSR  SIGMA_GAUSS'
;push,header,'#                   deg       deg      K       K       K       km/s      km/s'
;writeline,file,header,/prepend


;gass2.in on stream
;-------------------
restore,'gass2.dat'
add_tag,gstruc,'subcube',0L,gstruc
gstruc.subcube = 2
add_tag,gstruc,'X',0L,gstruc
add_tag,gstruc,'Y',0L,gstruc
gstruc.y = gstruc.lon+1
gstruc.x = gstruc.lat+2
newgstruc = gstruc
add_tag,gstruc,'MLON',0.0
add_tag,gstruc,'MLAT',0.0
gstruc.mlon = gstruc.glon
gstruc.mlat = gstruc.glat
mag2gal,gstruc.mlon,gstruc.mlat,glon,glat
gstruc.glon = glon
gstruc.glat = glat
; This cube spans -52.0 < MLON < -8.0
; Use -49.92 to -10.0
gd = where(gstruc.mlon gt -50.0 and gstruc.mlon le -10.0,ngd)
newgstruc=[newgstruc,gstruc[gd]]

;; Now write the ASCII file
;id = lindgen(n_elements(gstruc))+1
;fmt = '(I6,I5,I5,2F9.3,2F8.3,F9.3,F9.3,F9.3)' 
;;fmt = '(I8,2F7.1,2F8.3,F9.3,F9.3,F9.3)' 
;file = 'gass2.txt'
;writecol,file,id,gstruc.x,gstruc.y,gstruc.glon,gstruc.glat,gstruc.rms,gstruc.noise,$
;  gstruc.par[0],gstruc.par[1],gstruc.par[2],fmt=fmt
;
;undefine,header
;push,header,'#   NUM   X    Y   GLON      GLAT     RMS    NOISE   HEIGHT    VLSR  SIGMA_GAUSS'
;push,header,'#                   deg       deg      K       K       K       km/s      km/s'
;writeline,file,header,/prepend


;gass3.in on stream
;-------------------
restore,'gass3.dat'
add_tag,gstruc,'subcube',0L,gstruc
gstruc.subcube = 3
add_tag,gstruc,'X',0L,gstruc
add_tag,gstruc,'Y',0L,gstruc
gstruc.y = gstruc.lon+1
gstruc.x = gstruc.lat+2
newgstruc = gstruc
add_tag,gstruc,'MLON',0.0
add_tag,gstruc,'MLAT',0.0
gstruc.mlon = gstruc.glon
gstruc.mlat = gstruc.glat
mag2gal,gstruc.mlon,gstruc.mlat,glon,glat
gstruc.glon = glon
gstruc.glat = glat
; This cube spans -12.0 < MLON < 32.0
; Use -9.92 to +30.0
gd = where(gstruc.mlon gt -10.0 and gstruc.mlon le 30.0,ngd)
newgstruc=[newgstruc,gstruc[gd]]

;; Now write the ASCII file
;id = lindgen(n_elements(gstruc))+1
;fmt = '(I6,I5,I5,2F9.3,2F8.3,F9.3,F9.3,F9.3)' 
;;fmt = '(I8,2F7.1,2F8.3,F9.3,F9.3,F9.3)' 
;file = 'gass3.txt'
;writecol,file,id,gstruc.x,gstruc.y,gstruc.glon,gstruc.glat,gstruc.rms,gstruc.noise,$
;  gstruc.par[0],gstruc.par[1],gstruc.par[2],fmt=fmt
;
;undefine,header
;push,header,'#   NUM   X    Y   GLON      GLAT     RMS    NOISE   HEIGHT    VLSR  SIGMA_GAUSS'
;push,header,'#                   deg       deg      K       K       K       km/s      km/s'
;writeline,file,header,/prepend


;gass4.in on stream
;-------------------
restore,'gass4.dat'
add_tag,gstruc,'subcube',0L,gstruc
gstruc.subcube = 4
add_tag,gstruc,'X',0L,gstruc
add_tag,gstruc,'Y',0L,gstruc
gstruc.y = gstruc.lon+1
gstruc.x = gstruc.lat+2
newgstruc = gstruc
add_tag,gstruc,'MLON',0.0
add_tag,gstruc,'MLAT',0.0
gstruc.mlon = gstruc.glon
gstruc.mlat = gstruc.glat
mag2gal,gstruc.mlon,gstruc.mlat,glon,glat
gstruc.glon = glon
gstruc.glat = glat
; This cube spans 28.0 < MLON < 72.0
; Use +30.08 to 72.0
gd = where(gstruc.mlon gt 30.0,ngd)
newgstruc=[newgstruc,gstruc[gd]]

;; Now write the ASCII file
;id = lindgen(n_elements(gstruc))+1
;fmt = '(I6,I5,I5,2F9.3,2F8.3,F9.3,F9.3,F9.3)' 
;;fmt = '(I8,2F7.1,2F8.3,F9.3,F9.3,F9.3)' 
;file = 'gass4.txt'
;writecol,file,id,gstruc.x,gstruc.y,gstruc.glon,gstruc.glat,gstruc.rms,gstruc.noise,$
;  gstruc.par[0],gstruc.par[1],gstruc.par[2],fmt=fmt
;
;undefine,header
;push,header,'#   NUM   X    Y   GLON      GLAT     RMS    NOISE   HEIGHT    VLSR  SIGMA_GAUSS'
;push,header,'#                   deg       deg      K       K       K       km/s      km/s'
;writeline,file,header,/prepend

; Now combine them


;stop

gstruc = newgstruc

save,gstruc,file='gass.dat'


;; Now write the ASCII file
;id = lindgen(n_elements(gstruc))+1
;fmt2 = '(I6,I4,I5,I5,2F9.3,2F8.3,F9.3,F9.3,F9.3)' 
;file = 'gass.txt'
;writecol,file,id,gstruc.subcube,gstruc.x,gstruc.y,gstruc.glon,gstruc.glat,gstruc.rms,gstruc.noise,$
;  gstruc.par[0],gstruc.par[1],gstruc.par[2],fmt=fmt2
;undefine,header2
;push,header2,'#   NUM  SUB  X    Y   GLON      GLAT     RMS    NOISE   HEIGHT    VLSR  SIGMA_GAUSS'
;push,header2,'#                       deg       deg      K       K       K       km/s      km/s'
;writeline,file,header2,/prepend



stop

end
