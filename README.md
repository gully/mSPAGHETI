mSPAGHETI
=========
author: Michael Gully-Santiago (aka gully)
coauthor: Amanda Turbyfill
Date: January 2014
written in IDL

monochromatic Spectral Purity Analysis of Gratings with HDR to Enhance Technologies and  Instruments


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




