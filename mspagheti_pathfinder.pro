;monochromatic Spectral Purity Analysis of Gratings 
;    with HDR to Enhance Technologies & Instruments


;Function 0 Inputs
  ;Read in some user input for fits header
  ;History file?
  ;Mirror or grating?
function mspagheti_inputs, fn

dir1='/Users/gully/IDLWorkspace/mSPAGHETI/'
fn='mSPAGHETI_analysis_log_2014.csv'

;2) Read in the .xyz file data
tag_arr= ['VERSION','DATASTART','DELIMITER','MISSINGVALUE','COMMENTSYMBOL',$
   'FIELDCOUNT','FIELDTYPES','FIELDNAMES','FIELDLOCATIONS','FIELDGROUPS']
tag_fmt='F,J,B,F,A,J,J(18),A(18),J(18),J(18)'
create_struct, at_new, '',tag_arr, tag_fmt  

at_new.version=1.0
at_new.datastart=0
at_new.delimiter=44
at_new.missingvalue=!Values.F_NaN
at_new.COMMENTSYMBOL=';'
at_new.FIELDCOUNT=18
;print, strcompress((transpose(transpose(string(this_template.fieldtypes))+replicate(',', 18))))
at_new.FIELDTYPES=[ 3,  7,  7,  3,  7,  7,  3,  3,  4,  3,  4,  3,  7,  7,  7,  3,  3,  3]
at_new.FIELDNAMES=["Num", "Project", "Part",  "Sub_area",  "Directory", "Basename",  "Acq_Date",  $
      "d_mm", "wl",  "FL",  "ang", "N", "xc","yc","Note",  "f_missing", "f_nonstandard", "f_irror"]
at_new.FIELDLOCATIONS=[ 0,  3,  12,  15,  17,  76,  103,  112,  115,  121,  125,  137,  141,  142,  143,  168,  170,  172]
at_new.FIELDGROUPS=indgen(at_new.FIELDCOUNT)

d0=read_ascii(dir1+fn, template=at_new)

return, d0

print, 1

end


;Function 0.9, find (most of) the approximate x, y centers automatically
function mspagheti_axc_ayc, fn
d1=readfits(fn, /silent)
g=gauss2dfit(d1, A)
axc=a[4]
ayc=a[5]
return, [axc, ayc]
end

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
      
      mc_loopprogress,i,1,N-1
      
      ;Single frames
      fn=filename+string(i,FORMAT="(I03)")+'L.fit'
      if i gt 999 then fn =filename+string(i,FORMAT="(I04)")+'L.fit'
      wtemp=readfits(fn, /silent)
      darkname=filename+string(i,FORMAT="(I03)")+'D.fit'
      if i gt 999 then darkname=filename+string(i,FORMAT="(I04)")+'D.fit'
      darktemp=readfits(darkname, /silent)
      
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
 
      fn_out=f_basename+'xy.txt'
      forprint, cen_list[*, 0], cen_list[*, 1], TEXTOUT = fn_out, COMMENT = '; x, y (pixels)'
      return, fn_out

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
readcol, xy_fn, x_cen, y_cen

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

      mc_loopprogress,i,1,N-1      
      ;Single frames
      
      fn=f_basename+string(i,FORMAT="(I03)")+'L.fit'
      if i gt 999 then fn =f_basename+string(i,FORMAT="(I04)")+'L.fit'
      wtemp=readfits(fn, /silent)
      darkname=f_basename+string(i,FORMAT="(I03)")+'D.fit'
      if i gt 999 then darkname=f_basename+string(i,FORMAT="(I04)")+'D.fit'
      darktemp=readfits(darkname, /silent)
      
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
    fn_out=f_basename+'comb.fits'
    writefits,f_basename+'comb.fits',wcomb_n
    
    return, fn_out
    
end  
  
;------------------------------------
;Main procedure: mSPAGHETI_pathfinder
;------------------------------------
pro mSPAGHETI_pathfinder


;TODO
;We need a method to automatically assign the approximate x-y centers.
 
      d0=mspagheti_inputs(fn)


id=60

;Inputs:
dir1=d0.directory[id]
cd, dir1
filename=d0.basename[id]
d_mm=d0.d_mm[id]
fl_mm=d0.fl[id]
wl_nm=d0.wl[id]
bl_ang_deg=d0.ang[id]
N_frames=d0.N[id]
axc=float(d0.xc[id])
ayc=float(d0.yc[id])


;Function 1 Find drift 
      fn_xy=mspagheti_drift(filename, N_frames, axc, ayc)

;Function 2 Coadd files
      fn_comb=mspagheti_coadd(filename, fn_xy)
      ;Put in a flag for non-monotonic SNR     


;Todo:
;  Put some outputs into the log file:
  ;Combined file name
  ;achieved SNR
  ;achieved FWHM
  ;final xc yc of comb file
  ;

print, 1

end
