pro combine_gauss

; gass1
;------
undefine,gass1

; there's always overlap in ONE LAT.

;; gass1-1
;restore,'gass1-1.dat'
;add_tag,gstruc,'grid','gass1-1',gstruc
;push,gass1,gstruc
;
;; gass1-2
;restore,'gass1-2.dat'
;add_tag,gstruc,'grid','gass1-2',gstruc
;gd = where(gstruc.lat gt 0,ngd)
;push,gass1,gstruc[gd]
;	
;; gass1-3 
;restore,'gass1-3.dat'
;add_tag,gstruc,'grid','gass1-3',gstruc
;gd = where(gstruc.lat gt 0,ngd)
;push,gass1,gstruc[gd]
;
;; gass1-4 
;restore,'gass1-4.dat'
;add_tag,gstruc,'grid','gass1-4',gstruc
;gd = where(gstruc.lat gt 0,ngd)
;push,gass1,gstruc[gd]
;
;; gass1-5 
;restore,'gass1-5.dat'
;add_tag,gstruc,'grid','gass1-5',gstruc
;gd = where(gstruc.lat gt 0,ngd)
;push,gass1,gstruc[gd]

;save,gass1,file='gass1.dat'

;stop

; gass2
;------
;undefine,gass2
;
;; gass2-1
;restore,'gass2-1.dat'
;add_tag,gstruc,'grid','gass2-1',gstruc
;push,gass2,gstruc
;
;; gass2-2
;restore,'gass2-2.dat'
;add_tag,gstruc,'grid','gass2-2',gstruc
;gd = where(gstruc.lat gt 0,ngd)
;push,gass2,gstruc[gd]
;	
;; gass2-3 
;restore,'gass2-3.dat'
;add_tag,gstruc,'grid','gass2-3',gstruc
;gd = where(gstruc.lat gt 0,ngd)
;push,gass2,gstruc[gd]
;
;; gass2-4 
;restore,'gass2-4.dat'
;add_tag,gstruc,'grid','gass2-4',gstruc
;gd = where(gstruc.lat gt 0,ngd)
;push,gass2,gstruc[gd]
;
;; gass2-5 
;restore,'gass2-5.dat'
;add_tag,gstruc,'grid','gass2-5',gstruc
;gd = where(gstruc.lat ge 60,ngd)
;push,gass2,gstruc[gd]
;
;; gass2-5b
;restore,'gass2-5b.dat'
;add_tag,gstruc,'grid','gass2-5b',gstruc
;gd = where(gstruc.lat gt 0 and gstruc.lat lt 60,ngd)
;push,gass2,gstruc[gd]
;
; NEED TO ADD IN gass2-5b STIL!!!

;save,gass2,file='gass2.dat'

;stop

; gass3
;------
undefine,gass3

;gass3-1.dat	gass3-1c.dat	gass3-1e.dat	gass3-2b.dat	gass3-3.dat	gass3-5.dat
;gass3-1b.dat	gass3-1d.dat	gass3-2.dat	gass3-2c.dat	gass3-4.dat


;; gass3-1, 0-27
;restore,'gass3-1.dat'
;add_tag,gstruc,'grid','gass3-1',gstruc
;gd = where(gstruc.lat le 20,ngd)
;push,gass3,gstruc[gd]
;
;; gass3-1b, 15-40
;restore,'gass3-1b.dat'
;add_tag,gstruc,'grid','gass3-1b',gstruc
;gd = where(gstruc.lat gt 20 and gstruc.lat le 39,ngd)
;push,gass3,gstruc[gd]
;
;; gass3-1c, 38-63
;restore,'gass3-1c.dat'
;add_tag,gstruc,'grid','gass3-1c',gstruc
;gd = where(gstruc.lat gt 39 and gstruc.lat le 62,ngd)
;push,gass3,gstruc[gd]
;
;; gass3-1d, 61-86
;restore,'gass3-1d.dat'
;add_tag,gstruc,'grid','gass3-1d',gstruc
;gd = where(gstruc.lat gt 62 and gstruc.lat le 85,ngd)
;push,gass3,gstruc[gd]
;
;; gass3-1e, 84-110
;restore,'gass3-1e.dat'
;add_tag,gstruc,'grid','gass3-1e',gstruc
;gd = where(gstruc.lat gt 85,ngd)
;push,gass3,gstruc[gd]
;
;; gass3-2
;restore,'gass3-2.dat'
;add_tag,gstruc,'grid','gass3-2',gstruc
;gd = where(gstruc.lat gt 0,ngd)
;push,gass3,gstruc[gd]
;
;; don't need gass3-2b, gass3-2c
;
;; gass3-3 
;restore,'gass3-3.dat'
;add_tag,gstruc,'grid','gass3-3',gstruc
;gd = where(gstruc.lat gt 0,ngd)
;push,gass3,gstruc[gd]
;
;; gass3-4 
;restore,'gass3-4.dat'
;add_tag,gstruc,'grid','gass3-4',gstruc
;gd = where(gstruc.lat gt 0,ngd)
;push,gass3,gstruc[gd]
;
;; gass3-5b, 0-61
;restore,'gass3-5b.dat'
;add_tag,gstruc,'grid','gass3-5b',gstruc
;gd = where(gstruc.lat gt 0 and gstruc.lat le 60,ngd)
;push,gass3,gstruc[gd]
;
;; gass3-5, 59-110
;restore,'gass3-5.dat'
;add_tag,gstruc,'grid','gass3-5',gstruc
;gd = where(gstruc.lat gt 60,ngd)
;push,gass3,gstruc[gd]
;
;; 0-58 missing from gass3-5, running gass3-5b now
;
;save,gass3,file='gass3.dat'

;stop

; gass4
;------
undefine,gass4

; gass4-1.dat	gass4-2.dat	gass4-3.dat	gass4-4.dat	gass4-4b.dat	gass4-5.dat	gass4-5b.dat

;; gass4-1
;restore,'gass4-1.dat'
;add_tag,gstruc,'grid','gass4-1',gstruc
;push,gass4,gstruc
;
;; gass4-2
;restore,'gass4-2.dat'
;add_tag,gstruc,'grid','gass4-2',gstruc
;gd = where(gstruc.lat gt 0,ngd)
;push,gass4,gstruc[gd]
;	
;; gass4-3 
;restore,'gass4-3.dat'
;add_tag,gstruc,'grid','gass4-3',gstruc
;gd = where(gstruc.lat gt 0,ngd)
;push,gass4,gstruc[gd]
;
;; gass4-4 
;restore,'gass4-4.dat'
;add_tag,gstruc,'grid','gass4-4',gstruc
;gd = where(gstruc.lat gt 0,ngd)
;push,gass4,gstruc[gd]
;
;; don't need gass4-4b
;
;; gass4-5c, 0-6
;restore,'gass4-5c.dat'
;add_tag,gstruc,'grid','gass4-5c',gstruc
;gd = where(gstruc.lat gt 0 and gstruc.lat lt 6,ngd)
;push,gass4,gstruc[gd]
;
;; gass4-5d, 5-11
;restore,'gass4-5d.dat'
;add_tag,gstruc,'grid','gass4-5d',gstruc
;gd = where(gstruc.lat ge 6 and gstruc.lat lt 11,ngd)
;push,gass4,gstruc[gd]
;
;; gass4-5e, 10-16
;restore,'gass4-5e.dat'
;add_tag,gstruc,'grid','gass4-5e',gstruc
;gd = where(gstruc.lat ge 11 and gstruc.lat lt 16,ngd)
;push,gass4,gstruc[gd]
;
;; gass4-5f, 15-21
;restore,'gass4-5f.dat'
;add_tag,gstruc,'grid','gass4-5f',gstruc
;gd = where(gstruc.lat ge 16 and gstruc.lat lt 21,ngd)
;push,gass4,gstruc[gd]
;
;; gass4-5g, 20-26
;restore,'gass4-5g.dat'
;add_tag,gstruc,'grid','gass4-5g',gstruc
;gd = where(gstruc.lat ge 21 and gstruc.lat lt 26,ngd)
;push,gass4,gstruc[gd]
;
;; gass4-5h, 25-31
;restore,'gass4-5h.dat'
;add_tag,gstruc,'grid','gass4-5h',gstruc
;gd = where(gstruc.lat ge 26 and gstruc.lat lt 31,ngd)
;push,gass4,gstruc[gd]
;
;; gass4-5i, 30-36
;restore,'gass4-5i.dat'
;add_tag,gstruc,'grid','gass4-5i',gstruc
;gd = where(gstruc.lat ge 31 and gstruc.lat lt 36,ngd)
;push,gass4,gstruc[gd]
;
;; gass4-5j, 35-41
;restore,'gass4-5j.dat'
;add_tag,gstruc,'grid','gass4-5j',gstruc
;gd = where(gstruc.lat ge 36 and gstruc.lat lt 41,ngd)
;push,gass4,gstruc[gd]
;
;; gass4-5k, 40-46
;restore,'gass4-5k.dat'
;add_tag,gstruc,'grid','gass4-5k',gstruc
;gd = where(gstruc.lat ge 41 and gstruc.lat lt 46,ngd)
;push,gass4,gstruc[gd]
;
;; gass4-5l, 45-51
;restore,'gass4-5l.dat'
;add_tag,gstruc,'grid','gass4-5l',gstruc
;gd = where(gstruc.lat ge 46 and gstruc.lat lt 51,ngd)
;push,gass4,gstruc[gd]
;
;; gass4-5m, 50-56
;restore,'gass4-5m.dat'
;add_tag,gstruc,'grid','gass4-5m',gstruc
;gd = where(gstruc.lat ge 51 and gstruc.lat lt 56,ngd)
;push,gass4,gstruc[gd]
;
;; gass4-5n, 55-60
;restore,'gass4-5n.dat'
;add_tag,gstruc,'grid','gass4-5n',gstruc
;gd = where(gstruc.lat ge 56 and gstruc.lat lt 60,ngd)
;push,gass4,gstruc[gd]
;
;; gass4-5, 59-93
;restore,'gass4-5.dat'
;add_tag,gstruc,'grid','gass4-5',gstruc
;gd = where(gstruc.lat ge 60 and gstruc.lat le 92,ngd)
;push,gass4,gstruc[gd]
;
;; gass4-5b, 90-110
;restore,'gass4-5b.dat'
;add_tag,gstruc,'grid','gass4-5b',gstruc
;gd = where(gstruc.lat gt 92,ngd)
;push,gass4,gstruc[gd]
;
;save,gass4,file='gass4.dat'


; Combine all four pieces
;------------------------

restore,'gass1.dat'
gd = where(gass1.glon lt -50,ngd)
gass = gass1[gd]
undefine,gass1
restore,'gass2.dat'
gd = where(gass2.glon ge -50 and gass2.glon lt -10,ngd)
push,gass,gass2[gd]
undefine,gass2
restore,'gass3.dat'
gd = where(gass3.glon ge -10 and gass3.glon lt 30.0,ngd)
push,gass,gass3[gd]
undefine,gass3[gd]
restore,'gass4.dat
gd = where(gass4.glon ge 30.0,ngd)
push,gass,gass4[gd]
undefine,gass4

;save,gass,file='gass.dat'

; Need to use 0.0799d for the GLON/GLAT bin/step size

stop

end
