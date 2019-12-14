pro gass_plots

restore,'gass.dat'

gd = where(gass.glon gt -40 and gass.glon lt 15 and gass.glat gt -35 and gass.glat lt 10,ngd)
str = gass[gd]
undefine,gass
area = garea(str.par)

hess,str.glon,str.par[1],area,dx=0.0799,dy=1,/xflip,/log,yr=[10,400],top=50
; the MSLON/MSLAT step is 0.0799

ghess,str,'glon','cen','glat',dx=0.0799,dy=2,/xflip,xr=[-40,5],yr=[20,400],/total,top=2000,cut='glat gt -6'

g = where(str.glon le 5 and str.glat gt -6 and str.glat lt 10)
hess,str[g].glon,str[g].par[1],area[g],im,dx=0.0799,dy=2,yr=[20,400],/total,top=2000,xarr=xarr,yarr=yarr

!p.font = 0
setdisp
file = 'gassga_filaments'
ps_open,file,/color,thick=4,/encap
posim=[0.08,0.08,0.97,0.87]
poscol=[0.08,0.97,0.97,0.99]
displayc,im,xarr,yarr,xtit='L!dMS (deg)',ytit='V!dLSR!n (km/s)',tit='Total Intensity of GASS Gaussians (-6<B!dMS!n<10 deg)',/xflip,$
         charsize=1.2,posim=posim,poscol=poscol
; arrows for LMC filament
arrow,-8.2005,322.53993,-7.4580966,238.00623,/data,hthick=2.0,thick=4.0,/solid,co=255
arrow,-16.8,320,-15.773,228,/data,hthick=2.0,thick=4.0,/solid,co=255
arrow,-21.95-1,180.335+80,-21.95,180.335,/data,hthick=2.0,thick=4.0,/solid,co=255
arrow,-30.35-1,195.83+80,-30.35,195.83,/data,hthick=2.0,thick=4.0,/solid,co=255
arrow,-36.47-1,124.796+80,-36.47,124.796,/data,hthick=2.0,thick=4.0,/solid,co=255

plot,[0],/nodata,xr=reverse(minmax(xarr)),yr=minmax(yarr),xs=1,ys=1,/noerase,co=255,xtit=' ',ytit=' ',xtickformat='(A1)',ytickformat='(A1)',$
     charsize=1.2,position=posim
ps_close
ps2jpg,file+'.eps',/eps


stop

end
