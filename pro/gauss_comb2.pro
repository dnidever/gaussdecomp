pro gauss_comb2

; This program combines all hi7?.dat and hi8?.dat gstruc structures with
; the previous ones and transforms the coordinates to my "Magellanic"
; coordinates

; MAKE SURE TO GET RID OF THE **OVERLAP** !!!!!!
; There is a 4 degree overlap b/w all of strips.  Maybe
; I should start/stop them in the middle of the overlap
; 2 degrees to the left one and 2 degrees to the right one.
;;

;hi7a.in on rigel     lonr=[240,270], latr=[-20,40]
;hi7b.in on ultrafly  lonr=[266,296], latr=[-20,40]
;hi7c.in on helios    lonr=[292,322], latr=[-20,40]
;
;hi8a.in on calypso  [240,250], [-5,40], /NOBACK -
;hi8b.in on calypso  [248,258], [-5,40], /NOBACK -
;hi8c.in on beavis   [256,266], [-5,40], /NOBACK -
;hi8d.in on solstice [264,274], [-5,40], /NOBACK -
;hi8e.in on pulsar   [272,282], [-5,40], /NOBACK -
;hi8f.in on pleiades [280,290], [-5,40], /NOBACK -
;hi8g.in on pulsar   [288,298], [-5,40], /NOBACK -
;hi8h.in on mimas    [296,306], [-5,40], /NOBACK -
;hi8i.in on tethys   [304,315], [-5,40], /NOBACK -


;hi7a.in on rigel     lonr=[240,270], latr=[-20,40]
restore,'hi7a.dat'
gd1 = where(gstruc.lat gt -20.0 and gstruc.lat lt -5.0 and gstruc.lon ge 240.0 $
            and gstruc.lon le 268.0, ngd1)
newgstruc = gstruc(gd1)

;stop

;hi7b.in on ultrafly  lonr=[266,296], latr=[-20,40]
restore,'hi7b.dat'
gd2 = where(gstruc.lat gt -20.0 and gstruc.lat lt -5.0 and gstruc.lon gt 268.0 $
            and gstruc.lon le 294.0, ngd2)
newgstruc=[newgstruc,gstruc(gd2)]

;stop

;hi7c.in on helios    lonr=[292,322], latr=[-20,40]
restore,'hi7c.dat'
gd3 = where(gstruc.lat gt -20.0 and gstruc.lat lt -5.0 and gstruc.lon gt 294.0 $
            and gstruc.lon le 322.0, ngd2)
newgstruc=[newgstruc,gstruc(gd3)]

;stop

;hi8a.in on calypso  [240,250], [-5,40], /NOBACK -
restore,'hi8a.dat'
gd4 = where(gstruc.lat ge -5.0 and gstruc.lat le 40.0 and gstruc.lon ge 240.0 $
            and gstruc.lon le 249.0, ngd2)
newgstruc=[newgstruc,gstruc(gd4)]

;stop

;hi8b.in on calypso  [248,258], [-5,40], /NOBACK -
restore,'hi8b.dat'
gd5 = where(gstruc.lat ge -5.0 and gstruc.lat le 40.0 and gstruc.lon gt 249.0 $
            and gstruc.lon le 257.0, ngd2)
newgstruc=[newgstruc,gstruc(gd5)]

;stop

;hi8c.in on beavis   [256,266], [-5,40], /NOBACK -
restore,'hi8c.dat'
gd6 = where(gstruc.lat ge -5.0 and gstruc.lat le 40.0 and gstruc.lon gt 257.0 $
            and gstruc.lon le 265.0, ngd2)
newgstruc=[newgstruc,gstruc(gd6)]

;stop

;hi8d.in on solstice [264,274], [-5,40], /NOBACK -
restore,'hi8d.dat'
gd7 = where(gstruc.lat ge -5.0 and gstruc.lat le 40.0 and gstruc.lon gt 265.0 $
            and gstruc.lon le 273.0, ngd2)
newgstruc=[newgstruc,gstruc(gd7)]

;stop

;hi8e.in on pulsar   [272,282], [-5,40], /NOBACK -
restore,'hi8e.dat'
gd8 = where(gstruc.lat ge -5.0 and gstruc.lat le 40.0 and gstruc.lon gt 273.0 $
            and gstruc.lon le 281.0, ngd2)
newgstruc=[newgstruc,gstruc(gd8)]

;stop

;hi8f.in on pleiades [280,290], [-5,40], /NOBACK -
restore,'hi8f.dat'
gd9 = where(gstruc.lat ge -5.0 and gstruc.lat le 40.0 and gstruc.lon gt 281.0 $
            and gstruc.lon le 289.0, ngd2)
newgstruc=[newgstruc,gstruc(gd9)]

;stop

;hi8g.in on pulsar   [288,298], [-5,40], /NOBACK -
restore,'hi8g.dat'
gd10 = where(gstruc.lat ge -5.0 and gstruc.lat le 40.0 and gstruc.lon gt 289.0 $
            and gstruc.lon le 297.0, ngd2)
newgstruc=[newgstruc,gstruc(gd10)]

;stop

;hi8h.in on mimas    [296,306], [-5,40], /NOBACK -
restore,'hi8h.dat'
gd11 = where(gstruc.lat ge -5.0 and gstruc.lat le 40.0 and gstruc.lon gt 297.0 $
            and gstruc.lon le 305.0, ngd2)
newgstruc=[newgstruc,gstruc(gd11)]

;stop

;hi8i.in on tethys   [304,315], [-5,40], /NOBACK -
restore,'hi8i.dat'
gd12 = where(gstruc.lat ge -5.0 and gstruc.lat le 40.0 and gstruc.lon gt 305.0 $
            and gstruc.lon le 315.0, ngd2)
newgstruc=[newgstruc,gstruc(gd12)]


gstruc = newgstruc

save,gstruc,file='hi8.dat.orig'

lon = gstruc.lon
lat = gstruc.lat

gal2mag,lon,lat,mlon,mlat
;rotate_lb,lon,lat,[8.5,7.5],[285.5,-32.5],nlon,nlat,/noprint
;rotate_lb,lon,lat,[8.2679,7.5305],[285.6977,-44.3683],nlon,nlat,/noprint

; Making a new structure, adding the magellanic coordinates
dum = {mlon:0.,mlat:0.,glon:0.,glat:0.,rms:0.,noise:0.,par:fltarr(3),sigpar:fltarr(3)}
ng = n_elements(gstruc)
newgstruc = replicate(dum,ng)

newgstruc.glon = gstruc.lon
newgstruc.glat = gstruc.lat
newgstruc.rms = gstruc.rms
newgstruc.noise = gstruc.noise
newgstruc.par = gstruc.par
newgstruc.sigpar = gstruc.sigpar

newgstruc.mlon = mlon
newgstruc.mlat = mlat

gstruc = newgstruc

; Wrapping the longitudes to continuous negative values
bd = where(gstruc.mlon gt 200.,nbd)
gstruc(bd).mlon = gstruc(bd).mlon-360.

save,gstruc,file='hi8.dat'

stop

end
