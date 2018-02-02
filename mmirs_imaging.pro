;;
;; Hack to reduce imaging data by cutting and pasting
;; parts of  I. Chiligarian's MMIRS reduction
;; pipeline version 1.1 of 2015-01-22
;;
function get_filename,str,suffix=suffix
   if(n_elements(suffix) ne 1) then suffix=''
   return,[strsplit(str,',',/extract)+suffix]
end


pro mmirs_imaging,logfile, linear=linear,keepfirst=keepfirst,$
    verbose=verbose,debug=debug,biasframe=biasframe,badamp=badamp,$
    crosstalk=crosstalk,compress=compress,tmpdir=tmpdir,clean=clean

suffix='.fix.fits'

log = readlog(logfile)

rawdir = sxpar(log,'RAW_DIR')
rdir = sxpar(log,'R_DIR')
wdir = sxpar(log,'W_DIR')
;file_mkdir,wdir
spawn,'mkdir -p '+rdir
spawn,'mkdir -p '+wdir

grism = sxpar(log,'GRISM')
filter = sxpar(log,'FILTER')
slit_id = def_slit(logfile)

rawext = sxpar(log,'RAWEXT',count=cntre)
if(cntre eq 1) then suffix=suffix+rawext else rawext=''

sci_pref = sxpar(log,'SCI')
dithpos = sxpar(log,'DITHPOS')
dark_sci_all = sxpar(log,'DARKSCI')

sci_pref2 = sxpar(log,'SCI2',count=cntval)
diffmode = (cntval eq 1)? 1 : 0

arc_all = sxpar(log,'ARC')
dark_arc_all = sxpar(log,'DARKARC')

flat_all = sxpar(log,'FLAT')
dark_flat_all = sxpar(log,'DARKFLAT')

proc_stages=bytarr(9)
for i=0,n_elements(proc_stages)-1 do proc_stages[i]=sxpar(log,'S'+string(i+1,format='(i2.2)')+'PROC')

;; save directories and different images and image types associated
;; with these exposure in a structure

if(proc_stages[0] eq 1) then begin 
   conf_str = {rawdir:rawdir,$
               rdir:rdir,$
               wdir:wdir,$
               sci:get_filename(sci_pref),$
               sci_dark:get_filename(dark_sci_all),$
               arc:get_filename(arc_all),$
               arc_dark:get_filename(dark_arc_all),$
               flat:get_filename(flat_all),$
               flat_dark:get_filename(dark_flat_all),$
               misc:[''],misc_dark:get_filename(dark_sci_all)}
   
   if(diffmode eq 1) then begin
      conf_str.misc=[sci_pref2]
   endif
;;
;; The pipeline code seems to use the first dark frame as a bias
;; for the linearity correction; would it not be better to use
;; the median dark ?
;;
   biasframe = rawdir+conf_str.sci_dark[0]+'.fits'
;;
;; Verify if the linearised data cube exists, otherwise perform the
;; linearity correction
;;
;; This 
;; 1. creates a bias image from first 2 readouts of dark
;;
;; For each readout:
;;    2. subtract the bias 
;;    3. calculate and subtract the resistant mean of science image
;;       reference pixels 
;;    4. correct for detector cross-talk (if keyword is set)
;;    5. write results into image cube
;;
;; 6. correct non-linearity
;; 7. remove cosmic rays (if keyword "clean" is set)
;;    and flag saturated pixels
;; 

;;
;; this needs fixing in main code by ICh., removing the "+rawext"
;; otherwise the linearity correction is re-calculated for every
;; science image.
;;
   if(~file_test(conf_str.rdir+conf_str.sci+suffix)) then begin
      print,'processing ', conf_str.sci+suffix
      mmfixen_nonlin,badamp=badamp,conf_str.rawdir+conf_str.sci+'.fits',conf_str.rdir+conf_str.sci+suffix,biasframe=biasframe,compress=compress,verbose=verbose,crosstalk=crosstalk
;;   if(~file_test(conf_str.rdir+conf_str.sci+suffix+rawext)) then
;;   mmfixen_nonlin,badamp=badamp,conf_str.rawdir+conf_str.sci+'.fits',conf_str.rdir+conf_str.sci+suffix,biasframe=biasframe,compress=compress,verbose=verbose,crosstalk=crosstalk
   endif
;;
;; Same for the individual darks
;;
   for i=0,n_elements(conf_str.sci_dark)-1 do begin
      if(~file_test(conf_str.rdir+conf_str.sci_dark[i]+suffix)) then mmfixen_nonlin,badamp=badamp,conf_str.rawdir+conf_str.sci_dark[i]+'.fits',conf_str.rdir+conf_str.sci_dark[i]+suffix,biasframe=conf_str.rawdir+conf_str.sci_dark[i]+'.fits',compress=compress,verbose=verbose,crosstalk=crosstalk
;;      if(~file_test(conf_str.rdir+conf_str.sci_dark[i]+suffix+rawext)) then mmfixen_nonlin,badamp=badamp,conf_str.rawdir+conf_str.sci_dark[i]+'.fits',conf_str.rdir+conf_str.sci_dark[i]+suffix,biasframe=conf_str.rawdir+conf_str.sci_dark[i]+'.fits',compress=compress,verbose=verbose,crosstalk=crosstalk
   endfor
endif

print,'Subtracting dark frames'
;;
;; this part was lifted from subtract_dark_mmirs.pro
;; In the case of imaging, this will be the last step
;; and the output file name should not be hard-coded
;;
version=get_mmirs_pipeline_version()

if(n_elements(suffix) ne 1) then suffix='.fix.fits'

if((strlen(conf_str.sci[0]) gt 0) and (strlen(conf_str.sci_dark[0]) gt 0)) then begin
   h_sci=headfits(conf_str.rdir+conf_str.sci+suffix)
;;    sxaddpar,h_sci,'SOFTWARE',version,' data reduction software'
;;    if(n_elements(h_sci) eq 1) then message,'FITS file: '+conf_str.rdir+conf_str.sci+suffix+' not found'
   print, conf_str.rdir+conf_str.sci+suffix
    sci=mrdfits(conf_str.rdir+conf_str.sci+suffix,1,h_sci_ext,/silent)
;;    print,h_sci_ext
;;    sxaddpar,h_sci_ext,'BZERO',0.0
;;    sxaddpar,h_sci_ext,'SOFTWARE',version,' data reduction software'
    sxaddpar,h_sci_ext,'BZERO',0.0
    sxaddpar,h_sci_ext,'SOFTWARE',version,' data reduction software'
    dark=average_img(conf_str.rdir+conf_str.sci_dark+suffix,median=median)
    writefits,conf_str.wdir+conf_str.sci[0]+'_obj_dark.fits',0,h_sci
    print,conf_str.wdir+conf_str.sci[0]+'_obj_dark.fits'
    mwrfits,sci-dark,conf_str.wdir+conf_str.sci[0]+'_obj_dark.fits',h_sci_ext
endif

print,'done - status OK'

end
