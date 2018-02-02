pro correct_telluric,logfile,telluric_numb,diffmode=diffmode,noextr=noextr

log=readlog(logfile)
wdir=def_wdir(logfile)
n_tel=n_elements(telluric_numb)

for i=0,n_tel-1 do begin
    tellspec=mrdfits(wdir+'correction_tel_'+string(telluric_numb[i],format='(i2.2)')+'.fits',1,h,/silent)
    atm_model_tab=mrdfits(wdir+'correction_tel_'+string(telluric_numb[i],format='(i2.2)')+'.fits',2)
    if(i eq 0) then begin
        tell_arr=dblarr(n_elements(tellspec),n_tel)
        airmass_arr=dblarr(n_tel)
        water_vapor_arr=dblarr(n_tel)
    endif
    tell_arr[*,i]=tellspec
    airmass_arr[i]=sxpar(h,'AIRMASS')
    water_vapor_arr[i]=sxpar(h,'WVAPOR_M')
endfor

mask=read_mask_mmirs_ms(wdir+'mask_mos.txt')
n_slit=n_elements(mask)

obj_types=['obj','obj-sky']
if(keyword_set(diffmode)) then obj_types=[obj_types,'obj_diff','obj_diff-sky']

h=headfits(wdir+obj_types[i]+'_slits_lin.fits')
targ_airmass=sxpar(h,'AIRMASS')
targ_water_vapor=water_vapor_arr[0]

message,/inf,'Correcting telluric absorption for AIRMASS='+string(targ_airmass,format='(f7.3)')
tell_int=interp_telluric_correction(tell_arr,airmass_arr,targ_airmass,$
    water_vapor_arr,targ_water_vapor,atm_model_tab=atm_model_tab)

for i=0,n_elements(obj_types)-1 do begin
    print,'Correcting '+obj_types[i]+' for telluric absorption'
    file_in_2d=wdir+obj_types[i]+'_slits_lin.fits'
    file_out_2d=wdir+obj_types[i]+'_slits_corr.fits'

    if(not(keyword_set(noextr))) then begin
;;        file_in_1d=wdir+obj_types[i]+'_slits_extlr.fits'
;; cnaw: typo above ? 2017-07-20
        file_in_1d=wdir+obj_types[i]+'_slits_extr.fits'
        file_out_1d=wdir+obj_types[i]+'_slits_extr_corr.fits'
    endif

    hpri=headfits(file_in_2d)
    sxaddhist,'Corrected for telluric absorption',hpri
    writefits,file_out_2d,0,hpri

    for s=0,n_slit-1 do begin
        data_slit=mrdfits(file_in_2d,s+1,hext,/silent)
        ny_cur=sxpar(hext,'NAXIS2')
        for y=0,ny_cur-1 do data_slit[*,y]=data_slit[*,y]/tell_int
        mwrfits,data_slit,file_out_2d,hext
    endfor

    if(not(keyword_set(noextr))) then begin
       print,'correct_telluric: file_in_1d: ',file_in_1d
       hpri=headfits(file_in_1d)
       sxaddhist,'Corrected for telluric absorption',hpri
       writefits,file_out_1d,0,hpri
       data_extr=mrdfits(file_in_1d,1,hext,/silent)
       ny_cur=sxpar(hext,'NAXIS2')
       for y=0,ny_cur-1 do data_extr[*,y]=data_extr[*,y]/tell_int
       mwrfits,data_extr,file_out_1d,hext
    endif
endfor

end
