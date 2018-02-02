function get_nreads_ramp,inpfile,compress=compress,$
    exptime=exptime,ramp=ramp,header=h_img,gain=gain,rdnoise=rdnoise,pri_hdu=pri_hdu

h_pri=headfits(inpfile,compress=compress)
naxis=sxpar(h_pri,'NAXIS')
pri_hdu=0
if(naxis eq 2) then begin
    h_pri[0]='XTENSION= ''IMAGE''  /   /Primary FITS image array'
    sxaddpar,h_pri,'EXTNAME','IM'
    h_img=h_pri
    pri_hdu=1
endif else $
    h_img=headfits(inpfile,ext=1,compress=compress)
gain=sxpar(h_img,'GAIN')
rdnoise=sxpar(h_img,'RDNOISE')
nx=long(sxpar(h_img,'NAXIS1'))
ny=long(sxpar(h_img,'NAXIS2'))
exptime=sxpar(h_img,'EXPTIME')
sxdelpar,h_img,'BZERO'
sxaddpar,h_img,'BUNIT','counts/s'

if(pri_hdu eq 1) then begin
    ramp=(exptime eq 1)? 1 : 5
    n_img=1
endif else begin
    exptab=sxpar(h_img,'EXPTABLE')
    possec=strpos(exptab,'sec.tab')
    if(possec eq -1) then possec=strpos(exptab,'s.tab')
;; cnaw 2017-05-17 added this line :
    if(possec eq -1) then possec=strpos(exptab,'.tab')
    ramp=float(strmid(exptab,strpos(exptab,'ramp_')+5,$
        possec-strpos(exptab,'ramp_')-5))

    if(ramp gt exptime) then ramp=exptime
    if(ramp eq 4.5) then ramp=4.4919
    if(ramp eq 1.475) then ramp=1.72117
;    n_img=round(exptime/ramp)+1
    n_img=round(float(strmid(sxpar(h_img,'EXTNAME'),2)))
endelse

return,n_img

end
