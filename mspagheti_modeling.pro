;------prepare comparison data-------

;Function 4 Create model PSF

function mSPAGHETI_2Dairy, wl_nm, d_mm, px_scl

  k_um= 2.0*pi/(wl_nm/1000.0)

  arg= k_um*(d_mm*1000.0) ;k a
  
  ;make a radial circle, units of pixels
  xc = 512.0
  yc = 512.0
  rr=fltarr(xc*2, yc*2)
  for i = 0, xc*2-1 do begin
         for j = 0, yc*2-1 do begin
                 this_r=sqrt( (i-xc)^2 + (j-yc)^2 )
                 rr[i, j]=this_r
         endfor
  endfor
     
  ;convert radial circle to radian units:
  rr_rad=rr*px_scl
  
  full_arg=arg*sin(rr_rad)
   
  I_I0=(2.0*BeselJ(full_arg)/full_arg)^2.0
   
  ;Consider normalizing here...
   
  arr_out=I_I0
  return, arr_out

end


;Function 5 Determine if Mirror comparison exists
function mspagheti_mirexist, fn

;if file_test(fn) then begin
; arr=readfits(fn)
;else begin
; arr=fltarr(1024, 1024)+1.0E-8
;endif
;
;return, arr

return, 1

end


;Function 6 Determine if Zygo PSF exists
  ;Register Zygo PSF and Spectral Purity
function mspagheti_Zygexist, fn

;if file_test(fn) then do begin
; arr=readfits(fn)
;endif else begin
; arr=fltarr(1024, 1024)+1.0E-8
;endelse
;
;return, arr

print, 1

end  



;-------Analyze-------

;Function 7 2D intercomparisons saved to .eps files.
  ;Check for existing data (flags passed from above)
  
;function 7a- need a function to put all postage stamps on same spatial scale
function mspagheti_stamp_scl
print, 1

end
  
;function  7b
function mspagheti_2D_eps, dat_arr, mir_arr, zyg_arr, mod_arr, fn_out

device, decomposed=1
loadct, 13

SET_PLOT, 'PS'
!p.font=0
DEVICE, /ENCAPSULATED, /color, /helvetica, FILENAME = fn_out
device, xsize=12, ysize=12

!p.thick = 5
!x.thick = 5
!y.thick = 5
!z.thick = 5

     !P.Multi=[0,2,2]
     p = [0.02, 0.3, 0.98, 0.98]
     LoadCT, 13
     TVImage, bytscl(d1, -0.30, 0.30), Position=p, /keep_aspect_ratio, /nointerpolation
     Colorbar, Position=[p[0], p[1]-0.1, p[2], p[1]-0.05], range=[-0.3, 0.3], format='(f5.2)', divisions=2
     p = [0.02, 0.3, 0.98, 0.98]
     TVImage, bytscl(d2, -0.30, 0.30), Position=p, /keep_aspect_ratio, /nointerpolation
     Colorbar, Position=[p[0], p[1]-0.1, p[2], p[1]-0.05], range=[-0.3, 0.3], format='(f5.2)', divisions=2, ncolors=256
     p = [0.1, 0.3, 0.98, 0.98]
     LoadCT, 13
     TVImage, bytscl(d3, -0.1, 0.1), Position=p, /keep_aspect_ratio, /nointerpolation
     Colorbar, Position=[p[0], p[1]-0.1, p[2], p[1]-0.05], range=[-0.10, 0.10], format='(f5.2)', divisions=2
     !P.Multi =0




;plot, [0, 500], [0,500], thick=5, charthick=5, color='000000'xL, xrange=[0, 500],$
;yrange=[0, 500], /nodata, xtitle=xtit, ytitle=ytit, xstyle=1, ystyle=1
;tvimage, bytscl(d1, -0.2, 0.2), /keep_aspect_ratio, /nointerpolation
;contour, d,xv, yv, levels=[10, 20, 40, 80, 160, 320, 640], /overplot, /fill
;contour, d,xv, yv, levels=[13,14, 15, 16, 17, 20, 25, 30, 35, 40], /overplot, /fill

DEVICE, /CLOSE  ; Close the file.
SET_PLOT, 'X' ; Return plotting to Windows.
!p.font=-1




end  
  
 
;Function 8 1D cross-cut saved to .eps files.
  ;Check for existing data (flags passed from above)

;Function 9 Aperture photometry?


pro mSPAGHETI_modeling


print, 1


end