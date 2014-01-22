;monochromatic Spectral Purity Analysis of Gratings 
;    with HDR to Enhance Technologies & Instruments


;Function 0a Analysis Inputs
function mspagheti_inputs, fn

dir1='/Users/gully/IDLWorkspace/mSPAGHETI/'
fn='mSPAGHETI_analysis_log_2014.csv'

;at=ascii_template(dir1+fn)

;2) Read in the .xyz file data
tag_arr= ['VERSION','DATASTART','DELIMITER','MISSINGVALUE','COMMENTSYMBOL',$
   'FIELDCOUNT','FIELDTYPES','FIELDNAMES','FIELDLOCATIONS','FIELDGROUPS']
tag_fmt='F,J,B,F,A,J,J(19),A(19),J(19),J(19)'
create_struct, at_new, '',tag_arr, tag_fmt  

at_new.version=1.0
at_new.datastart=0
at_new.delimiter=44
at_new.missingvalue=!Values.F_NaN
at_new.COMMENTSYMBOL=';'
at_new.FIELDCOUNT=19
;print, strcompress((transpose(transpose(string(at.fieldtypes))+replicate(',', 19))))
at_new.FIELDTYPES=[ 3,  7,  7,  3,  7,  7,  3,  4,  4,  3,  4,  3,  4,  4,  7,  3,  3,  3,  3]
at_new.FIELDNAMES=["Num", "Project", "Part",  "Sub_area",  "Directory", "Basename",  "Acq_Date",  $
      "d_mm", "wl",  "FL",  "ang", "N", "xc","yc","Note",  "f_missing", "f_nonstandard", "f_mirror", "f_reduce"]
at_new.FIELDLOCATIONS=[0,  3,  12,  15,  17,  76,  103,  112,  117,  123,  127,  132,  136,  142,  148,  173,  175,  177,  179]
at_new.FIELDGROUPS=indgen(at_new.FIELDCOUNT)

d0=read_ascii(dir1+fn, template=at_new)

return, d0

end


;Function 0b Result Outputs
function mspagheti_outputs

dir1='/Users/gully/IDLWorkspace/mSPAGHETI/'
fn='mSPAGHETI_results_log_2014.csv'

;at=ascii_template(dir1+fn)

;2) Read in the .xyz file data
tag_arr= ['VERSION','DATASTART','DELIMITER','MISSINGVALUE','COMMENTSYMBOL',$
   'FIELDCOUNT','FIELDTYPES','FIELDNAMES','FIELDLOCATIONS','FIELDGROUPS']
tag_fmt='F,J,B,F,A,J,J(27),A(27),J(27),J(27)'
create_struct, at_new, '',tag_arr, tag_fmt  

at_new.version=1.0
at_new.datastart=9
at_new.delimiter=44
at_new.missingvalue=!Values.F_NaN
at_new.COMMENTSYMBOL=';'
at_new.FIELDCOUNT=27
;print, strcompress((transpose(transpose(string(at.fieldtypes))+replicate(',', 27))))
at_new.FIELDTYPES=[ 3,  3,  3,  3,  3,  3,  7,  4,  4,  4,  4,  4,  7,  7,  4,  4,  4,  4,  4,  7,  4,  3,  4,  4,  4,  4,  4]
at_new.FIELDNAMES=["Num", "range_x",  "range_y",  "stddev_x", "stddev_y", "drift_t",  "drift_fn", $
       "hist_peak",  "hist_sig", "log_SNR",  "SNR_slope",  "SNR_offset", "bad_frames", "coadd_fn",$
        "rot_ang",  "rot_ang_unc",  "xc_fin", "yc_fin", "rot_time", "rot_fn", "grat_period", $
         "obs_order",  "ang_adj_order",  "plt_scale",  "FWHM_pred",  "FWHM_meas",  "FWHM_lam_d"]
at_new.FIELDLOCATIONS=[ 0,  3,  6,  9,  12,  15,  18,  25,  30,  35,  40,  46,  52,  64,  71,  78,  84,  90,  98,  105,  113,  120,  123,  129,  136,  142,  148]
at_new.FIELDGROUPS=indgen(at_new.FIELDCOUNT)

d0=read_ascii(dir1+fn, template=at_new)

return, d0

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
      xx= cen_list[*, 0]
      yy=cen_list[*, 1]
      xrange=max(xx)-min(xx)
      yrange=max(yy)-min(yy)
      stddevx=stddev(xx)
      stddevy=stddev(yy)
      
      forprint, cen_list[*, 0], cen_list[*, 1], TEXTOUT = fn_out, COMMENT = '; x, y (pixels)'
      end_time = SYSTIME(1)
      tot_time = end_time - start_time ; in seconds 
      
      
      struct_out={drift_fn:fn_out, $
                  drift_t:tot_time, $
                  range_x:xrange,$
                  range_y:yrange,$
                  stddev_x:stddevx,$
                  stddev_y:stddevy}
      
      return, struct_out

end


;Function 1.9 Write to results log
  ;inputs: a structure of results from a process
  ;        the unique id number from the first column of the analysis and result logs
  ;        the structure of results
  ;outputs: updates the entries for that num_id
function mspagheti_populate_results, this_struct, num_id, res_struct

;First, find and list the tags in the result struct
;Second, get tags of structure
;Third, find the row where the num_id matches
;Fourth, find where the tags match the column names
;Fifth, populate the result_structure with values
;Sixth, save the database


;1) find and list the tags in the result struct
these_tags=get_tags(this_struct)

;2) get tags of structure
res_tags=get_tags(res_struct)

;3) find the row where the num_id matches
this_i=where(num_id eq res_struct.num)

;4) find where the tags match the column names
match, these_tags, res_tags, suba, subb

;5) populate the result_structure with values
n_matched=n_elements(suba)
if (n_matched ne -1) then begin
  for i=0, n_matched-1 do begin
    exec_string = 'res_struct'+these_tags[suba[i]]+'[num_id] = this_struct'+these_tags[suba[i]]
    temp1=execute(exec_string)
   endfor
endif else begin
  print, "The tags did not match"
endelse

return, res_struct

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
dark=wcomb
bf=4.0 ;binning factor
wcomb_e=fltarr(1024*bf, 1024*bf)


;Big for loop  
  N=n_elements(x_cen)
  SNR_arr=fltarr(N)
  bad_string=''
  
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
      wcomb=wcomb+w    
      ;estimate the noise properties
      histogauss,wcomb, a1, /nofit, charsize=3, /noplot
      SNR_arr[i]=max(wcomb)/a1[2]
      if SNR_arr[i] lt SNR_arr[i-1] then bad_string=strcompress(bad_string+'/'+string(i))
          
     ;Expand the data by rebinning by a binning factor, bf
      exp_dat=rebin(w, 1024*bf, 1024*bf) ;expanded data
      
      ;Shift the now expanded data by dx, dy
      sh_ex_dat=shift(exp_dat, x_cen[i]*bf, y_cen[i]*bf) ;shifted, expanded data
      
      ;big w combined
      wcomb_e=wcomb_e+sh_ex_dat
      
    endfor

    wcomb=rebin(wcomb_e, 1024, 1024)
    wcomb_n=wcomb/max(wcomb)
    histogauss,wcomb, a_out, /nofit, charsize=3
    
    ;Fit for the coefficient
    n1=alog10(1.0+findgen(N-1))
    y1=alog10(SNR_arr[1:*])
    coeff=linfit(n1, y1, /double)
    forprint, n1, y1, TEXTOUT = f_basename+'SNR.txt', COMMENT = '; N, SNR'
    
    ;Write the non-rotated frame
    fn_out=f_basename+'coadd.fits'
    writefits,fn_out,wcomb_n

    end_time=systime(1)
    tot_time=end_time-start_time
    
    struct_out={coadd_fn:fn_out,$
                hist_peak:a_out[1],$
                hist_sig:a_out[2],$
                log_SNR:alog10(1.0/a_out[2]),$
                SNR_slope:coeff[1],$
                SNR_offset:coeff[0],$
                bad_frames:bad_string}
    
    return, struct_out
    
end  
  
;------------------------------------
;Main procedure: mSPAGHETI_pathfinder
;------------------------------------
pro mSPAGHETI_pathfinder

 
d0=mspagheti_inputs(fn)

;Loop over the targets you wish to reduce
;The f_reduce flag should be set to 1 for these cases.

N_sources=n_elements(d0.f_reduce)
for i=0, N_sources-1 do begin

  id=i
  
  if d0.f_reduce[id] then begin
    ;Inputs:
    dir1=d0.directory[id]
    cd, dir1
    filename=d0.basename[id]
    d_mm=d0.d_mm[id]
    fl_mm=d0.fl[id]
    wl_nm=d0.wl[id]
    bl_ang_deg=d0.ang[id]
    N_frames=13;d0.N[id]
    axc=float(d0.xc[id])
    ayc=float(d0.yc[id])
      
      
    ;Function 1 Find drift 
          drift_out=mspagheti_drift(filename, N_frames, axc, ayc)
          
res_struct=mspagheti_outputs()
temp1=mspagheti_populate_results(drift_out, id, res_struct)
    
    ;Function 2 Coadd files
          struct_out=mspagheti_coadd(filename, drift_out.drift_fn)
          ;Put in a flag for non-monotonic SNR
    print, '-------------'
    print, 'Reducing '+strcompress(string(id))
    print, '-------------'
  endif else begin
    print, 'Skipping '+strcompress(string(id))
  endelse
endfor
  
    ;TODO
    ;look for OBOE in coadd and xy functions.
 
print, 1

end
