function generate_extr_krnl,dithpos,e_apw,e_meth=e_meth,$
    ccdscl=ccdscl,c_krnl=c_krnl,n_krnl=n_krnl

if(n_elements(e_meth) ne 1) then e_meth=1
if(n_elements(ccdscl) ne 1) then ccdscl=0.2

n_dith=n_elements(dithpos)

c_krnl = fix(max(abs(dithpos))/ccdscl) + fix(e_apw*3)
n_krnl = 2*c_krnl+1

extr_krnl = dblarr(n_krnl)

for i=0,n_dith-1 do begin
    if(e_meth eq 2) then begin
        d = psf_gaussian(npix=[n_krnl],ndim=1,$
            centroid=[c_krnl+dithpos[i]/ccdscl],fwhm=[e_apw],/nor,/double)
    endif else begin
        d = dblarr(n_krnl)
        d[c_krnl-floor(e_apw/2.0):c_krnl+floor(e_apw/2.0)]=1.0
        d0=d ;;;/total(d)
        d=shift_s(d0,dithpos[i]/ccdscl)
    endelse
    extr_krnl[*] = extr_krnl[*] + d*((-1.0)^i)
endfor

return,extr_krnl
end


pro extract_1d_slits,logfile,obstype,$
    dith_from_box=dith_from_box, box_exp=box_exp, $
    diffmode=diffmode, detect=detect, optimal=optimal, error=error,$
    suffix=suffix,writesuffix=writesuffix

log=readlog(logfile)
wdir=sxpar(log,'W_DIR')
if(n_elements(suffix) ne 1) then suffix='_lin'
if(n_elements(writesuffix) ne 1) then writesuffix=''

;; extraction method : (1) for boxcar; (2) for optimal extraction
e_meth = sxpar(log,'EXTMETH')
if(keyword_set(optimal)) then e_meth=2 ;;;; forcing optimal extraction

;; extraction aperture in pixels
e_apw = sxpar(log,'EXTAPW')

val=sxpar(log,'DITHPOS',count=cntval)
dithpos=(cntval eq 1)? double(strsplit(val,',',/extract)) : 0.0

if(keyword_set(diffmode)) then begin
    val2=sxpar(log,'DITHPOS2',count=cntval2)
    dithpos2=(cntval2 eq 1)? double(strsplit(val2,',',/extract)) : 0.0
    dithpos=[dithpos,dithpos2]
endif

ccdscl=0.2 ;;; 0.2 arcsec per pix -- hardcoded CCD scale

mask=read_mask_mmirs_ms(wdir+'mask_mos.txt')
n_slits=n_elements(mask)

inpfile=wdir+obstype+'_slits'+suffix+'.fits'
hpri=headfits(inpfile)

if(keyword_set(dith_from_box)) then begin
    if(n_elements(box_exp) ne 1) then box_exp='obj-sky'
    boxlist=where(mask.type eq 'BOX',cbox)
    for i=0,cbox-1 do begin
        t=mrdfits(wdir+box_exp+'_slits'+suffix+'.fits',boxlist[i]+1,hb,/silent)
        val=sxpar(hb,'DITHCOMP',count=cntval)
        dithpos_tmp=(cntval eq 1)? double(strsplit(val,',',/extract)) : 0.0
        if(i eq 0) then dithpos_arr=dblarr(n_elements(dithpos_tmp),cbox)
        dithpos_arr[*,i]=dithpos_tmp[*]/double(cbox)
    endfor
    if(keyword_set(diffmode)) then begin
        d_dithpos=dithpos[1]-dithpos[0]
    endif
    dithpos=(cbox gt 0)? total(dithpos_arr,2) : [0.0]
    if(keyword_set(diffmode)) then dithpos=[dithpos,dithpos+d_dithpos]
    print,'Using empirically estimated dithering positions: ',dithpos
endif

extr_krnl=generate_extr_krnl(dithpos,e_apw,ccdscl=ccdscl,e_meth=e_meth,$
    c_krnl=c_krnl,n_krnl=n_krnl)*e_apw

outfile=wdir+obstype+'_slits_extr'+writesuffix+'.fits'
writefits,outfile,0,hpri
for i=0,n_slits-1 do begin
    print,'Extracting from slit #',i+1,'/',n_slits
;; cnaw
    print,'extract_1d_slits: inpfile ', inpfile
    slit_img=mrdfits(inpfile,i+1,himg,/silent)
    nx=sxpar(himg,'NAXIS1')
    ny=sxpar(himg,'NAXIS2')
    if(i eq 0) then out_img=dblarr(nx,n_slits)

    if(keyword_set(detect)) then begin
       print,'extract_1d_slits: keyword_set(detect) is true'
        obj_pos=sxpar(himg,'OBJPOS',count=cnt_obj)
        obj_fwhm=sxpar(himg,'OBJFWHM')
        if(cnt_obj eq 1) then print,'Object extracted in slit '+string(i+1),'  pos=',obj_pos,' fwhm=',obj_fwhm else begin
            obj_pos=0.0
            obj_fwhm=e_apw
        endelse
        
        extr_krnl=(e_meth eq 2)? generate_extr_krnl([obj_pos],obj_fwhm,ccdscl=ccdscl,e_meth=e_meth,$
            c_krnl=c_krnl,n_krnl=n_krnl) : generate_extr_krnl([obj_pos],e_apw,ccdscl=ccdscl,e_meth=e_meth,$
            c_krnl=c_krnl,n_krnl=n_krnl)

;        endif else begin
;            extr_krnl=generate_extr_krnl([0.0],e_apw,ccdscl=ccdscl,e_meth=e_meth,$
;                c_krnl=c_krnl,n_krnl=n_krnl)
;        endelse
    endif
;; cnaw
    print,'extract_1d_slits: ny, n_krnl', ny, n_krnl
    krnl_cur=dblarr(ny)
    c_k_c=(ny-1)/2
;; cnaw
    print,'extract_1d_slits: c_k_c',c_k_c

    if(ny ge n_krnl) then begin
        krnl_cur[c_k_c-c_krnl:c_k_c+c_krnl]=extr_krnl
    endif else begin
        krnl_cur[*]=extr_krnl[c_krnl-c_k_c:c_krnl+c_k_c]
    endelse
    
;; cnaw
    print,'extract_1d_slits: n_elements(krnl_cur)',n_elements(krnl_cur),krnl_cur
    goodk = where(abs(krnl_cur) gt 1e-4)
;; cnaw
    print,'extract_1d_slits:n_elements(goodk), goodk',n_elements(goodk), goodk
    for x=0,nx-1 do begin
        out_img[x,i]=(keyword_set(error))? $
            sqrt(total((transpose(slit_img[x,goodk])*krnl_cur[goodk])^2,nan=diffmode)) : $
            total(transpose(slit_img[x,goodk])*krnl_cur[goodk],nan=diffmode)
    endfor

endfor
sxdelpar,himg,'YOFFSET'
sxaddpar,himg,'DITHPOS',sxpar(log,'DITHPOS')
if(n_slits eq 1) then sxdelpar,himg,'NAXIS2'

delkw=['']
for i=0,n_elements(himg)-1 do begin
    if(strmid(himg[i],0,4) eq 'MASK') then delkw=[delkw,strcompress(strmid(himg[i],0,8),/remove)]
endfor
n_bkw=n_elements(delkw)
if(n_bkw gt 1) then for i=1,n_bkw-1 do sxdelpar,himg,delkw[i]

mwrfits,out_img,outfile,himg
mwrfits,mask,outfile

end
