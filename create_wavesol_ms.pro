pro create_wavesol_ms,LOGFILE,slit=slit,ndeg=ndeg,smy=smy,$
    debug=debug,plot=plot,y_ndeg=y_ndeg,oh=oh

aaa=''
if(n_elements(ndeg) ne 1) then ndeg=3
if(n_elements(y_ndeg) ne 1) then y_ndeg=3
wdir=def_wdir(LOGFILE)

arc_file=(keyword_set(oh))? 'obj_ff.fits' : 'arc_ff.fits'
arc_image=readfits(wdir+arc_file,h)
dist_map=mrdfits(wdir+'dist_map.fits',1)

arc_r=poly_2d(rotate(arc_image,5),dist_map.kxwrap,dist_map.kywrap,1)
arc_image=rotate(arc_r,5)
arc_r=0

Nx=sxpar(h,'NAXIS1')
Ny=sxpar(h,'NAXIS2')
mask=read_mask_mmirs_ms(wdir+'mask_mos.txt')

y_scl = 30.2101

n_slits=n_elements(mask)
if(n_elements(slit) eq 0) then slit=findgen(n_slits)
n_slits=n_elements(slit)

y_slits=get_slit_region(mask,Nx=Nx,Ny=Ny,ratio=1.0,dist_map=dist_map,slit_geom=slit_geom)

pix_mask_data=replicate({y_pix:0,x_mask:!values.f_nan,y_mask:!values.f_nan,$
    w_pix:!values.f_nan,h_pix:!values.d_nan,slit:-1,$
    wl_sol:dblarr(Ndeg+1)+!values.d_nan,wl_s_err:!values.d_nan},Ny)
tags=tag_names(slit_geom[0])
for i=0,n_elements(tags)-1 do pix_mask_data.(i)=slit_geom.(i)
disper_all=dblarr(Ndeg+2,Ny)

;if(file_test(wdir+'disper.fits') and (not(keyword_set(overwrite)))) then begin
;    d=readfits(wdir+'disper.fits')
;    if(array_equal(size(disper_all),size(d))) then disper_all=d
;    d=0
;endif
;
;if(file_test(wdir+'disper_table.fits') and (not(keyword_set(overwrite)))) then begin
;    d=mrdfits(wdir+'disper_table.fits',1)
;    if(array_equal(size(pix_mask_data),size(d))) then begin
;        if(n_elements(pix_mask_data[0].wl_sol) eq n_elements(d[0].wl_sol)) then pix_mask_data=d
;    endif
;    d=0
;endif

for j=0,n_slits-1 do begin
    i=slit[j]
    ymin_cur = y_slits[0,i]
    ymax_cur = y_slits[1,i]

    if(mask[i].type ne 'TARGET') then continue

    print,'Processing slit: ',i+1,'/',n_slits,ymin_cur,ymax_cur
    fwhm=(mask[i].width*y_scl) > 2.0
    status_slit=crea_disper_mmirs_slit(logfile,i,fwhm,Ndeg,yfit='Yfit',plot=plot,NdegY=y_ndeg,$
        xpos=mask[i].x,ypos=mask[i].y,smy=smy,pssuff='_slit'+string(i+1,format='(i2.2)'),$
        ymin_out=ymin_cur,ymax_out=ymax_cur,oh=oh)
    if(status_slit eq 0) then begin
        dtmp=readfits(wdir+'disper.fits',/silent)
        disper_all[*,ymin_cur:ymax_cur]=dtmp[*,*]
        pix_mask_data[ymin_cur:ymax_cur].wl_sol=dtmp[0:Ndeg,*]
        pix_mask_data[ymin_cur:ymax_cur].wl_s_err=transpose(dtmp[Ndeg+1,*])
        if(keyword_set(debug)) then begin
            aaa=''
            print,'npix=',ymax_cur-ymin_cur+1,'/',$
                n_elements(where(finite(pix_mask_data[ymin_cur:ymax_cur].y_mask) eq 1))
            read,aaa
        endif
    endif
endfor
;b_d=where(finite(pix_mask_data.wl_s_err and pix_mask-data) ne 1)
;pix_mask_data[b_d].x_mask=!values.f_nan
;pix_mask_data[b_d].y_mask=!values.f_nan

writefits,wdir+'disper.fits',disper_all
writefits,wdir+'disper_table.fits',0
mwrfits,pix_mask_data,wdir+'disper_table.fits'

end
