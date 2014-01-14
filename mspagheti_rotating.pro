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
      TVImage, bytscl(disp_data, -5.0, -2.0), /nointerpolation
      
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
      TVImage, bytscl(disp_data, -5.0, -2.0), /nointerpolation
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
    fn_out=file_basename(fn_comb, '.fits')+'_rot.fits'
    writefits,fn_out,wcn_rot
    return, fn_out
    
end    




pro mSPAGHETI_rotating


;Function 3a Determine rotation
      theta_deg=mspagheti_detRot(fn_comb)

;Function 3b Optionally rotate
      fn_rot=mspagheti_perfRot(fn_comb, theta_deg)

end
