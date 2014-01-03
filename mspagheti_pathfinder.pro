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
      forprint, cen_list[*, 0], cen_list[*, 1], TEXTOUT = f_basename, FORMAT = format, SILENT = SILENT, $ 
      STARTLINE = startline, NUMLINE = numline, COMMENT = comment, $
      SUBSET = subset, NoCOMMENT=Nocomment,STDOUT=stdout
;+
end

;Function 2 Coadd files
  ;Read in files and x,y centes
  ;Coadd files
  ;Determine image properties reached (per pixel)
  ;Report image properties reached (per pixel)
  ;Make fits header
  ;Save combined fits file

;Function 3 Optionally rotate
  ;Determine rotation (by hand?)
  ;Perform rotation


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