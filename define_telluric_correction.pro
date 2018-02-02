function bin_hr_spec,wl_hr,spec_hr,wl_lr
    n_wl=n_elements(wl_lr)
    spec_lr=dblarr(n_wl)
    spec_hr_cnt=histogram(wl_hr,min=wl_lr[0]-(wl_lr[1]-wl_lr[0])/2.0,$
                                max=wl_lr[n_wl-1]+(wl_lr[n_wl-1]-wl_lr[n_wl-2])/2.0,$
                                binsize=(wl_lr[1]-wl_lr[0]),reverse_ind=ind_spec_hr)
    n_bins=(n_elements(spec_hr_cnt) < n_wl)

    for i=1L,n_bins-1L do $
        if((ind_spec_hr[i]-ind_spec_hr[i-1]) gt 0) then $
            spec_lr[i]=total(spec_hr[ind_spec_hr[ind_spec_hr[i-1]:ind_spec_hr[i]-1L]])

    spec_lr=spec_lr/double(spec_hr_cnt)

    return, spec_lr
end


function get_shift_ccorr,vec1,vec2,poly_cont=poly_cont,nmag=nmag,dpix=dpix,maxccorr=maxccorr,gauss=gauss
    n1=n_elements(vec1)
    n2=n_elements(vec2)
    maxccorr=!values.f_nan
    if(n1 ne n2) then begin
        message,/inf,'Vec1 and Vec2 have different lenghts. Returning NaN'
        return,!values.f_nan
    endif
    if(n_elements(nmag) ne 1) then nmag=10L
    if(n_elements(dpix) ne 1) then dpix=17

    x_cross=findgen(2*nmag*dpix+1)/double(nmag)-dpix

    x_orig=findgen(n1)
    vec1a=vec1
    vec2a=vec2
    bvec1a=where(finite(vec1a) ne 1,cbvec1a,compl=gvec1a,ncompl=cgvec1a)
    bvec2a=where(finite(vec2a) ne 1,cbvec2a,compl=gvec2a,ncompl=cgvec2a)
    if(cgvec1a lt 20 or cgvec2a lt 20) then return,!values.f_nan
    if(cbvec1a gt 0) then vec1a[bvec1a]=interpol(vec1a[gvec1a],x_orig[gvec1a],x_orig[bvec1a])
    if(cbvec2a gt 0) then vec2a[bvec2a]=interpol(vec2a[gvec2a],x_orig[gvec2a],x_orig[bvec2a])

    if(n_elements(poly_cont) eq 1) then begin
        kc1=poly_fit(x_orig,vec1a,poly_cont,yfit=c1)
        vec1a=vec1a-c1
        kc2=poly_fit(x_orig,vec2a,poly_cont,yfit=c2)
        vec2a=vec2a-c2
    endif

    vec1a_m = congrid(vec1a,n1*long(nmag),cubic=-0.5)
    vec2a_m = congrid(vec2a,n2*long(nmag),cubic=-0.5)
    c_cross=c_correlate(vec1a_m,vec2a_m,x_cross*double(nmag))
    maxccorr=max(c_cross,/nan,x_max)
    dx_max=x_cross[x_max]
    if(keyword_set(gauss)) then begin
       gau=gaussfit(x_cross,c_cross,G,nterms=3)
       dx_max=G[1]
    endif

    return,dx_max
end



pro define_telluric_correction,logfile,n_tel=n_tel,modelstar=modelstar,vr=vr,$
    absolute=absolute,Vmag_tell=Vmag_tell,debug=debug

log=readlog(logfile)
wdir=def_wdir(logfile)

if(n_elements(vr) ne 1) then vr=0.0
if(n_elements(modelstar) ne 1) then modelstar='a0v'

mask=read_mask_mmirs_ms(wdir+'mask_mos.txt')

n_slit=n_elements(mask)

if(n_elements(n_tel) ne 1) then n_tel=1

star_pref=(n_tel gt 0)? 'star_tel_'+string(n_tel,format='(i2.2)') : 'obj-sky'
file_2d=wdir+star_pref+'_slits_lin.fits'
file_flat=wdir+'flatn_slits_lin.fits'
file_1d=wdir+star_pref+'_slits_extr.fits'
file_out_corr=wdir+'correction_tel_'+string(n_tel,format='(i2.2)')+'.fits'
file_out_mean=wdir+'tel_'+string(n_tel,format='(i2.2)')+'_mean.fits'
file_out_norm=wdir+star_pref+'_norm.fits'

objpos_list=[-1]
active_slits=[-1]
inactive_slits=[-1]
obj_fwhm=[-1]
airmass=[-1]
for s=0,n_slit-1 do begin
    r=mrdfits(file_2d,s+1,h,/silent)
    objpos=sxpar(h,'OBJPOS',count=count)
    if(count eq 1) then begin
        objpos_list=[objpos_list,objpos]
        active_slits=[active_slits,s]
        obj_fwhm=[obj_fwhm,sxpar(h,'OBJFWHM')]
        airmass=[airmass,sxpar(h,'AIRMASS')]
    endif else inactive_slits=[inactive_slits,s]
endfor

n_active=n_elements(active_slits)-1
if(n_active eq 0) then begin
    message,/inf,'Telluric star not found in any of the slits. Cannot continue'
    return
endif

objpos_list=objpos_list[1:*]
active_slits=active_slits[1:*]
obj_fwhm=obj_fwhm[1:*]
airmass=airmass[1:*]

hpri=headfits(file_1d)
data_telluric=mrdfits(file_1d,1,h,/silent)
flatn_arr=data_telluric*0.0+1.0
parse_spechdr,h,wl=wl

if(keyword_set(absolute)) then begin
    for i=0,n_elements(active_slits)-1 do begin
        flatn=mrdfits(file_flat,active_slits[i]+1,/silent)
        flatn_arr[*,active_slits[i]]=flatn[*,n_elements(flatn[0,*])/2+objpos_list[i]]
    endfor
endif

n_wl=n_elements(wl)
data_telluric_orig=data_telluric
if(n_elements(inactive_slits) gt 1) then data_telluric_orig[*,inactive_slits[1:*]]=!values.f_nan

data_telluric=data_telluric[*,active_slits]
flatn_arr=flatn_arr[*,active_slits]

tel_norm_k=dblarr(n_active)


grism = def_grism(logfile,filter=filter)

if(grism eq 'H') then begin
    wl_min=1500.0
    wl_max=1700.0
    if(filter eq 'H') then begin
        wl_min_atm=1770.0
        wl_max_atm=1850.0
    endif else begin ;;; HK filter assumed
        wl_min_atm=1770.0
        wl_max_atm=1950.0
    endelse
endif

if(grism eq 'HK') then begin
    wl_min=1500.0
    wl_max=1700.0
    wl_min_atm = 1800.0
    wl_max_atm = 1950.0
endif

if(grism eq 'J') then begin
    wl_min=1180.0
    wl_max=1300.0
    wl_min_atm = 1150.0
    wl_max_atm = 1180.0
endif

if(grism eq 'H3000') then begin
    wl_min=1500.0
    wl_max=1750.0
    wl_min_atm=1755.0
    wl_max_atm=1782.0
endif

if(grism eq 'K3000') then begin
    wl_min=2000.0
    wl_max=2250.0
    wl_min_atm=1980.0
    wl_max_atm=2080.0
endif

tell_norm=data_telluric*0.0
norm_wl=where(wl ge wl_min and wl le wl_max, cnorm_wl)
absorb_wl=where(wl ge wl_min_atm and wl le wl_max_atm)

print,'Reading the atmosphere model'
;; cnaw 2017-07-20
calib_path = getenv('IDL_MMIRS')
;; cnaw 2017-07-20
atm_model=mrdfits(calib_path+'calib_MMIRS/sky_transmission/cptrans_zm_43_15.fits.gz',1,/silent)
atm_mmirs = bin_hr_spec(atm_model.wave,atm_model.transmission,wl)
wv_s=['23','43','76','100']
wv = [2.3,4.3,7.6,10.0]
am_s=['10','15','20']
am = [1.0,1.5,2.0]

atm_model_tab=replicate({water_vapor:0.0,airmass:0.0,wl:wl,transmission:atm_mmirs},12)
for i=0,3 do begin
    for j=0,2 do begin
;; cnaw 2017-07-20
        atm_model_cur=mrdfits(calib_path+'calib_MMIRS/sky_transmission/cptrans_zm_'+wv_s[i]+'_'+am_s[j]+'.fits.gz',1,/silent)
        atm_model_tab[j*4+i].water_vapor=wv[i]
        atm_model_tab[j*4+i].airmass=am[j]
        atm_model_tab[j*4+i].transmission=bin_hr_spec(atm_model_cur.wave,atm_model_cur.transmission,wl)
    endfor
endfor

print,'done'

wl_shift_pix=dblarr(n_active)

for i=0,n_active-1 do begin
    tel_norm_k[i]=(keyword_set(absolute))? 1d : median(data_telluric[norm_wl,i])
    tell_norm[*,i]=data_telluric[*,i]/tel_norm_k[i]
    if(keyword_set(absolute)) then tell_norm[*,i]=tell_norm[*,i]*flatn_arr[*,i]
    wl_shift_pix[i]=get_shift_ccorr(atm_mmirs[absorb_wl],tell_norm[absorb_wl,i],/gauss,maxcc=maxcc,poly_cont=4)
    print,'Shift for slit',active_slits[i]+1,' is ',wl_shift_pix[i],' pixels, max(corr)=',maxcc
    
    tell_norm[*,i]=shift_s(tell_norm[*,i],-wl_shift_pix[i])
    data_telluric_orig[*,active_slits[i]]=tell_norm[*,i]
endfor

writefits,file_out_norm,0,hpri
mwrfits,data_telluric_orig,file_out_norm,h

tell_med=(n_active gt 1)? median(tell_norm,dim=2) : tell_norm

print,'Reading the stellar model spectrum for '+modelstar
;; cnaw 2017-07-20
star_model=mrdfits(calib_path+'calib_MMIRS/telluric/'+modelstar+'.fits.gz',1,/silent)
;;;star_model.wave=star_model.wave/(1.0+vr/299792.5d)
print,'done'

star_mmirs=bin_hr_spec(star_model.wave,star_model.flux,wl)
if(keyword_set(absolute) and n_elements(Vmag_tell) eq 1) then begin
;; cnaw 2017-07-20
    v_b90 = read_asc(calib_path+'calib_MMIRS/telluric/V_B90.dat',1) ;;; Bessel 1990 V filter
    model_v=bin_hr_spec(star_model.wave,star_model.flux,transpose(v_b90[0,*])/10.0)
    model_flux_v=total(model_v*v_b90[1,*])/total(v_b90[1,*])
    k_norm = model_flux_v/3.63d-9*10d^(0.4*Vmag_tell)
    star_mmirs=star_mmirs/k_norm ;;; now the star is in F_lambda erg/cm^2/s/A
    mirror_area = 3.318e+5*0.89 ;;; 6.5-m mirror area in cm^2 minus the central obscuration (11%)
    star_mmirs=star_mmirs*mirror_area*5.034e+7*wl*10.0 ;;; erg/s/A to photons/s/A conversion, wavelength is in nm
    gain = 1.0 ;;; data already in e-/sec
    tell_med=tell_med*gain/(wl[1]-wl[0])/10.0 ;; wavelength is in nm
endif else begin
    star_mmirs=star_mmirs/median(star_mmirs[norm_wl])
endelse

wl_shift_vr = get_shift_ccorr(star_mmirs[norm_wl],tell_med[norm_wl],/gauss,poly_cont=4,maxcc=maxcc)
vr = ((wl_shift_vr*[wl[1]-wl[0]])/wl[cnorm_wl/2.0])*299792.458d
print,'Radial velocity shift:',wl_shift_vr,' pix = ',vr,' km/s, max(corr)=',maxcc

star_mmirs=shift_s(star_mmirs,wl_shift_vr)

;;;;; still to convolve with the instrumental response (not a simple gaussian)

fwhm_slits=mean(mask[active_slits].width)*30.2101
fwhm_meas =median(obj_fwhm)/0.2012

if(fwhm_slits gt fwhm_meas) then begin
    if(keyword_set(debug)) then print,'FWHM_slits>FWHM_meas',fwhm_slits,fwhm_meas
    krnl=psf_gaussian(ndim=1,npix=[15],fwhm=[fwhm_meas],/norm)
    corrfunc = tell_med/convol(star_mmirs,krnl)
    krnl_slit = dblarr(round(fwhm_slits))+1.0
    krnl_slit=krnl_slit/total(krnl_slit)
    corrfunc = convol(corrfunc,krnl_slit)
    for i=0,11 do begin
        t_vec=convol(atm_model_tab[i].transmission,krnl_slit)
        atm_model_tab[i].transmission=t_vec
    endfor
endif else begin
    if(keyword_set(debug)) then print,'FWHM_slits<FWHM_meas',fwhm_slits,fwhm_meas
    krnl0=psf_gaussian(ndim=1,npix=[15],fwhm=[2.0],/norm) ;;;; 2.0 pixel Gaussian
    krnl_slit = dblarr(round(fwhm_slits))+1.0
    krnl=convol(krnl0,krnl_slit/total(krnl_slit))
    krnl=krnl/total(krnl)
    corrfunc = tell_med/convol(star_mmirs,krnl)
    for i=0,11 do begin
        t_vec=convol(atm_model_tab[i].transmission,krnl)
        atm_model_tab[i].transmission=t_vec
    endfor
endelse

if(keyword_set(debug)) then print,krnl
;;;print,'Gaussian FWHM: ',fwhm,' pix'
airmass_cur=((mean(airmass)) > 1.0)
sxaddpar,hpri,'AIRMASS',airmass_cur,' mean airmass'

writefits,file_out_mean,0,hpri
sxdelpar,h,'NAXIS2'
mwrfits,tell_med,file_out_mean,h

writefits,file_out_corr,0,hpri
sxdelpar,h,'NAXIS2'
;;;sxaddpar,h,'MODFWHM',fwhm,'Kernel FWHM (pix) used to convolve the model star'

debug_val=(keyword_set(debug))? 1 : 0
atm_model_int = interp_telluric_model(atm_model_tab,airmass=airmass_cur)
water_vapor_bfit=tnmin('compare_atm_transmission',[3.0],$
       functargs={data:corrfunc,atm_model_tab:atm_model_int,$
                  grism:grism,filter:filter,debug:debug_val},$
       parinfo=[{limited:[1,1],limits:[min(atm_model_tab.water_vapor),max(atm_model_tab.water_vapor)]}],$
       /autoder)

water_vapor_value = ((water_vapor_bfit ge min(atm_model_tab.water_vapor)) and $
                     (water_vapor_bfit le max(atm_model_tab.water_vapor)))? $
        water_vapor_bfit[0] : 4.0

sxaddpar,h,'WVAPOR_M',water_vapor_value,'Model based estimated water vapor, mm'

mwrfits,corrfunc,file_out_corr,h


mwrfits,atm_model_tab,file_out_corr

end
