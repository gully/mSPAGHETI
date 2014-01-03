;monochromatic Spectral Purity Analysis of Gratings 
;    with HDR to Enhance Technologies & Instruments


;Function 0 Inputs
  ;Read in some user input for fits header
  ;History file?
  ;Mirror or grating?


;Function 1 Find drift 
  ;Align/register (Gaussian with sub frame and approx x,y-center)
  ;[experimental] Align/register (cross-correlation)
  ;track goodness of fit
        ;TODO: clean-up windows
        ;identify potential errors
  ;output: 
  ;  text file with xy centers
function mspagheti_drift, f_basename, N_frames, axc, ayc
  
;time the running time
start_time = SYSTIME(1) 

pi=3.141592654

filename=f_basename
f1=f_basename+string(N_frames/2,FORMAT="(I03)")+'L.fit'

;Combined files and single frames
  wcomb=fltarr(1024,1024)
  w=wcomb
  dark=wcomb
  bf=4.0 ;binning factor
  wcomb_e=fltarr(1024*bf, 1024*bf)
  

;Big for loop  
  N=N_frames
  cen_list=dblarr(N+1,2)
  par_lst=dblarr(N+1, 4)
  snr_trac=fltarr(N+1)
  
  for i=1,N-1 do begin
      
      ;Single frames
      fn=filename+string(i,FORMAT="(I03)")+'L.fit'
      if i gt 999 then fn =filename+string(i,FORMAT="(I04)")+'L.fit'
      wtemp=readfits(fn)
      darkname=filename+string(i,FORMAT="(I03)")+'D.fit'
      if i gt 999 then darkname=filename+string(i,FORMAT="(I04)")+'D.fit'
      darktemp=readfits(darkname)
      
      ;Subtract the dark
      w=float(wtemp)-float(darktemp)
      
      ;Search in a small box around the approx center
      ;fit a gaussian there to get precise center info
      p=40  
      t_sa=w[axc-p:axc+p, ayc-p:ayc+p] ;80 x 80 array for the x,y center fit
      d_sa=w[*, ayc-p:ayc+p-1] ; 1024 x 20 array for the rebinning and summing
      t_sa=t_sa-median(w)
      t_sa=t_sa/max(t_sa)

      ;Fit a 2d Gaussian
      difs=OB_fit_airy(t_sa, p+1, p+1, 2.4)
      
      cen_list[i-1,0]=difs[0]; max_x-axc
      cen_list[i-1, 1]=difs[1] ;max_y-ayc
      
      endfor
      forprint, cen_list[*, 0], cen_list[*, 1], TEXTOUT = f_basename, COMMENT = '; x, y (pixels)'
      return, 1

end

;Function 2 Coadd files
  ;Read in files and x,y centes
  ;Coadd files
  ;Determine image properties reached (per pixel)
  ;Report image properties reached (per pixel)
  ;Make fits header
  ;Save combined fits file
  
function mspagheti_coadd, f_basename, xy_fn

;Read in files and x,y centers
readcol, x_cen, y_cen, filename=xy_fn

  ;Coadd files
  ;Determine image properties reached (per pixel)
  ;Report image properties reached (per pixel)
  ;Make fits header
  ;Save combined fits file

;time the running time
start_time = SYSTIME(1) 

pi=3.141592654

;Combined files and single frames
wcomb=fltarr(1024,1024)
w=wcomb
dark=wcomb
bf=4.0 ;binning factor
wcomb_e=fltarr(1024*bf, 1024*bf)
  

;Big for loop  
  N=n_elements(x_cen)
  
  for i=1,N-1 do begin
      
      ;Single frames
      
      fn=f_basename+string(i,FORMAT="(I03)")+'L.fit'
      if i gt 999 then fn =f_basename+string(i,FORMAT="(I04)")+'L.fit'
      wtemp=readfits(fn)
      darkname=f_basename+string(i,FORMAT="(I03)")+'D.fit'
      if i gt 999 then darkname=f_basename+string(i,FORMAT="(I04)")+'D.fit'
      darktemp=readfits(darkname)
      
      ;Subtract the dark
      w=float(wtemp)-float(darktemp)
          
     ;Expand the data by rebinning by a binning factor, bf
      exp_dat=rebin(w, 1024*bf, 1024*bf) ;expanded data
      
      ;Shift the now expanded data by dx, dy
      sh_ex_dat=shift(exp_dat, x_cen[i]*bf, y_cen[i]*bf) ;shifted, expanded data
      
      ;big w combined
      wcomb_e=wcomb_e+sh_ex_dat
      
    endfor

    wcomb=rebin(wcomb_e, 1024, 1024)
    wcomb_n=wcomb/max(wcomb)
    
    ;Write the non-rotated frame
     writefits,f_basename+'comb.fits',wcomb_n
    
end  
  

;Function 3a Determine rotation
  ;Determine rotation by hand
  
function mspagheti_detRot, fn_comb  

      wcomb_n=readfits(fn_comb)

      ;finds position of maximum pixel
      maxval=max(wcomb_n, pos)
      xcen=pos mod 1024
      ycen=pos/1024
      
      device, decomposed=0
      loadct, 13
      window, 10, xsize=1024, ysize=1024, Title='Click on the left side order'
      disp_data=alog10(wcomb_n)
      TVImage, bytscl(disp_data, -4.0, 0.0), /nointerpolation
      
      ;selects a distant point on the spec. dimension
      ;  to the left of the brightest peak for angle determination
      wset, 10
      cursor, x1, y1, /down, /device
      xyouts, x1, y1, "X", /device, charsize=2
      
      wait, 0.4
      
      ;selects a distant point on the spec. dimension
      ;  to the right of the brightest peak for angle determination
      window, 10, xsize=1024, ysize=1024, Title='Click on the right side order'
      disp_data=alog10(wcomb_n)
      TVImage, bytscl(disp_data, -4.0, 0.0), /nointerpolation
      cursor, x2, y2, /down, /device
      xyouts, x2, y2, "X", /device, charsize=2
      
      
      ;Refine the position of the ghosts for rotation
      s=12
      ;Left side ghost
      gh1_sa=wcomb_n[x1-s:x1+s, y1-s:y1+s] ;10 x 10 array for the x,y center fit
      gh1_d=OB_fit_airy(gh1_sa, s, s, 1.65)
      
      ;Right side ghost
      gh2_sa=wcomb_n[x2-s:x2+s, y2-s:y2+s] ;10 x 10 array for the x,y center fit
      gh2_d=OB_fit_airy(gh2_sa, s, s, 1.65)
      
      ;Calculate the slope, m
      ;Fit to the PSF makes a refined estimate for the slope
      x1f=double(x1-gh1_d[0])
      y1f=double(y1-gh1_d[1])
      x2f=double(x2-gh2_d[0])
      y2f=double(y2-gh2_d[1])
      
      ;Slope, m
      m=double(y2f-y1f)/double(x2f-x1f)
      
      ;output the info
      print, 'The slope is:'
      print, m
      
      theta_rad=atan(m)
      theta_deg=180.0/3.141592654*theta_rad
      
      print, "theta (deg) is:", theta_deg
      
      return, theta_deg
end    

;Function 3b Optionally rotate
  ;Perform rotation
  
function mspagheti_perfRot, fn_comb, theta_deg

      wcomb_n=readfits(fn_comb)
      ;perform the rotation.  do not use /pivot)
      wcomb_ex=fltarr(1024*3.0, 1024*3.0)
      wcomb_ex[1024:1024+1024-1, 1024:1024+1024-1]=wcomb_n
      maxval=max(wcomb_ex, pos)
      xcen=pos mod (1024*3)
      ycen=pos/(1024*3)
      wcn_rot=rot(wcomb_ex, theta_deg, 1, xcen, ycen, /cubic)

    ;Write the fits files
    writefits,fn_comb+'_rot.fits',wcn_rot
    return, 1
    
end    


;------prepare comparison data-------

;Function 4 Create model PSF

;Function 5 Determine if Mirror comparison exists

;Function 6 Determine if Zygo PSF exists
  ;Register Zygo PSF and Spectral Purity

;-------Analyze-------

;Function 7 2D intercomparisons saved to .eps files.
  ;Check for existing data (flags passed from above)  

;Function 8 1D cross-cut saved to .eps files.
  ;Check for existing data (flags passed from above)

;Function 9 Aperture photometry?



pro mSPAGHETI_pathfinder

;Inputs:
dir1='/Users/gully/Astronomy/silicon/APRA_JPL/TJ07/D/'
cd, dir1
filename='TJ07_d6mm_f838mm_632nm-'
d_mm=6.0
fl_mm=838.0
wl_nm=632.8
bl_ang_deg=54.7

;Function 1 Find drift 
  ;Align/register (Gaussian with sub frame and approx x,y-center)
  ;inputs: 
      ;file basename
      ;Number of frames
      ;aproximate x center
      ;approximate y center
  ;outputs: 
      ;text file with xy centers vs frame number
      drift_out=mspagheti_drift(f_basename, N_frames, axc, ayc)


print, 1

end