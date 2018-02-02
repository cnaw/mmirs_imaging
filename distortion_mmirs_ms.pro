pro distortion_mmirs_ms,logfile,image_type,gzip=gzip,deg1=deg1,deg2=deg2,diff=diff,force=force

wdir=def_wdir(logfile)
grism=def_grism(logfile,filter=filter)
slit=def_slit(logfile)


if(n_params() eq 1) then image_type='flat'
suff=(keyword_set(gzip))? '.gz' : ''

period=get_mmirs_period(wdir+image_type+'_dark.fits'+suff)

if(slit ne 'mos') then begin ;;;; longslit assumed
    file_copy,'calib_MMIRS/LS/'+period+'/dist_map_'+grism+'_'+filter+'.fits',wdir+'/dist_map.fits',/overwrite
    return
endif

lin_params=linearisation_params(grism,filter)
if((grism eq 'H' and filter eq 'H') or (grism eq 'J' and filter eq 'J')) then diff=1

shortspec = ((grism eq 'H' and filter eq 'H') or (grism eq 'J' and filter eq 'J'))? 1 : 0

if(n_elements(deg1) ne 1) then deg1=(shortspec eq 1)? 2 : 3
if(n_elements(deg2) ne 1) then deg2=(shortspec eq 1)? 1 : 3

h=headfits(wdir+image_type+'_dark.fits'+suff)
img=transpose(mrdfits(wdir+image_type+'_dark.fits'+suff,1,himg,/silent))
Nx=sxpar(himg,'NAXIS2')
Ny=sxpar(himg,'NAXIS1')
;mask_data=mrdfits(wdir+'mask_def.fits',1)
mask=read_mask_mmirs_ms(wdir+'mask_mos.txt',logfile=logfile)
n_slit=n_elements(mask)

n_st=(shortspec eq 1)? 8L : 32L
wx=(shortspec eq 1)? 120 : 20
p_ref=median(img[Nx/2-wx:Nx/2+wx,*],dim=1)

skippix=0
nsegments=deg2 + 3
overlap=0.2
z=(keyword_set(diff))? 20.0 : 10.0
dy_amp = (shortspec eq 1)? 15.0 : 35.0
dy=-dy_amp+findgen(2.0*dy_amp*z+1L)/z

dyarr=dblarr(n_st,nsegments)
c_arr=dyarr
xgrid=dblarr(n_st)
ygrid=dblarr(nsegments)
for xx=fix(n_st/2),n_st-1 do begin
    xcur=(xx+0.5)*(double(Ny)/n_st)
    xgrid[xx]=xcur
    p_cur=median(img[xcur-wx:xcur+wx,*],dim=1)
    if(xx eq fix(n_st/2)) then begin
        xref_min=xcur-wx
        xref_max=xcur+wx
        prof_ref = p_cur
    endif
;    plot,p_cur,xs=1,title='x='+string(xcur) & oplot,p_ref,col=128 & zzz='' & read,zzz
;    plot,dy,dy*0,xs=1,ys=1,yr=[0,1],title='x='+string(xcur)
    if (keyword_set(diff) and (xx eq fix(n_st/2))) then begin
        p_ref = p_cur
        continue
    endif
    for i=0,nsegments-1 do begin
        ymin = (i eq 0)? skippix : (double(i)-overlap)*Ny/nsegments
        ymax = (i eq nsegments-1)? Ny-skippix-1 : (double(i)+1.0+overlap)*Ny/nsegments
        if(ymin lt 0) then ymin=0
        if(ymax ge Ny) then ymax=Ny-1
        ymin=fix(ymin)
        ymax=fix(ymax)
        ygrid[i]=(ymin+ymax)/2.0
        vec1=rebin(p_ref[ymin:ymax],(ymax-ymin+1L)*z)
        vec2=rebin(p_cur[ymin:ymax],(ymax-ymin+1L)*z)
        cf=c_correlate(vec1-min(vec1),vec2-min(vec2),dy*z)
;        oplot,dy,cf,col=30*(i+1)
        nn=max(cf,cnn)
        if((cnn eq 0) or (cnn eq (n_elements(dy)-1L)) or (nn lt 0.2)) then begin
            dyarr[xx,i]=!values.f_nan
            c_arr[xx,i]=!values.f_nan
        endif else begin
            dyarr[xx,i]=(keyword_set(diff))? dy[cnn] : dy[cnn]
            c_arr[xx,i]=nn
        endelse
    endfor

    if(keyword_set(diff)) then p_ref=p_cur
;    read,zzz
endfor

for xx=fix(n_st/2),0,-1 do begin
    xcur=(xx+0.5)*(double(Ny)/n_st)
    xgrid[xx]=xcur
    p_cur=median(img[xcur-wx:xcur+wx,*],dim=1)
;    if(xx eq fix(n_st/2)) then begin
;        prof_ref = p_cur
;    endif
;    plot,p_cur,xs=1,title='x='+string(xcur) & oplot,p_ref,col=128 & zzz='' & read,zzz
;    plot,dy,dy*0,xs=1,ys=1,yr=[0,1],title='x='+string(xcur)
    if (keyword_set(diff) and (xx eq fix(n_st/2))) then begin
        p_ref = p_cur
        continue
    endif
    for i=0,nsegments-1 do begin
        ymin = (i eq 0)? skippix : (double(i)-overlap)*Ny/nsegments
        ymax = (i eq nsegments-1)? Ny-skippix-1 : (double(i)+1.0+overlap)*Ny/nsegments
        if(ymin lt 0) then ymin=0
        if(ymax ge Ny) then ymax=Ny-1
        ymin=fix(ymin)
        ymax=fix(ymax)
        ygrid[i]=(ymin+ymax)/2.0
        vec1=rebin(p_ref[ymin:ymax],(ymax-ymin+1L)*z)
        vec2=rebin(p_cur[ymin:ymax],(ymax-ymin+1L)*z)
        cf=c_correlate(vec1,vec2,dy*z)
;        oplot,dy,cf,col=30*(i+1)
        nn=max(cf,cnn)
        if((cnn eq 0) or (cnn eq (n_elements(dy)-1L)) or (nn lt 0.2)) then begin
            dyarr[xx,i]=!values.f_nan
            c_arr[xx,i]=!values.f_nan
        endif else begin
            dyarr[xx,i]=(keyword_set(diff))? -dy[cnn] : dy[cnn]
            c_arr[xx,i]=nn
        endelse
    endfor

    if(keyword_set(diff)) then p_ref=p_cur
;    read,zzz
endfor

if(keyword_set(diff)) then begin
    dyarr[fix(n_st/2),*]=(dyarr[fix(n_st/2)+1,*]+dyarr[fix(n_st/2)-1,*])/2.0
    dyarr_orig=dyarr
    for i=0,nsegments-1 do begin
        gg=where(finite(dyarr_orig[*,i]) eq 1,cgg)
        k_t = robust_poly_fit(xgrid[gg],dyarr_orig[gg,i]/(xgrid[1]-xgrid[0]),2)
        k_t_int = [0,k_t]/[1.0,1.0,2.0,3.0]
        dyarr[*,i]=poly(xgrid,k_t_int)
    endfor
    dyarr=dyarr-median(dyarr[fix(n_st/2)-1:fix(n_st/2),*])
endif
;stop

xgrid2d=(xgrid # (dblarr(n_elements(ygrid))+1.0))
ygrid2d=((dblarr(n_elements(xgrid))+1.0) # ygrid)
g=where(finite(dyarr) eq 1,cg)
darr=dblarr(3,cg)
darr[0,*]=xgrid2d[g]
darr[1,*]=ygrid2d[g]
darr[2,*]=dyarr[g]

e=sfit_2deg(darr,deg1,deg2,kx=kx,/irr,/max)

polywarp,xgrid2d,ygrid2d,$
         xgrid2d,ygrid2d-poly2d(xgrid,ygrid,kx,deg1=deg1,deg2=deg2),$
         3,Kxwr,Kywr

polywarp,xgrid2d,ygrid2d-poly2d(xgrid,ygrid,kx,deg1=deg1,deg2=deg2),$
         xgrid2d,ygrid2d,$
         3,Kxwr_inv,Kywr_inv

;; added full path 
res_str_stored=mrdfits('~/idl/mmirs/pipeline/calib_MMIRS/MOS/'+period+'/dist_map_'+grism+'_'+filter+'_mos.fits',1,/silent)

res_str={XGRID:xgrid, YGRID:ygrid, $
         KXWRAP:kxwr,KYWRAP:kywr, $
         KXWRAP_INV:kxwr_inv,KYWRAP_INV:kywr_inv, $
         dyarr:dyarr, corrarr:c_arr,$
         kx_2dfit:kx,deg1:deg1,deg2:deg2,mask_dy0:-20.0,mask_y_scl:30.2101d}

if(grism eq 'HK' and filter eq 'HK' and n_elements(res_str_stored) eq 1) then begin
    u=[0,1,4,5,7,8,9,10,11,12,13,14,15]
    c2kywrap_inv=total(((res_str.kywrap_inv/res_str_stored.kywrap_inv)[u]-1)^2)
    if(c2kywrap_inv gt 2.0 and (not keyword_set(force))) then begin
        message,/inf,'Problem computing distorion map. Using stored values'
        res_str=res_str_stored
    endif
endif

if(grism eq 'J' and filter eq 'zJ' and n_elements(res_str_stored) eq 1) then begin
    u=[0,1,4,5,7,8,9,10,11,12,13,14,15]
    res_str_stored.mask_dy0=-20.0
    c2kywrap_inv=total(((res_str.kywrap_inv/res_str_stored.kywrap_inv)[u]-1)^2)
    if(c2kywrap_inv gt 2.0 and (not keyword_set(force))) then begin
        message,/inf,'Problem computing distorion map. Using stored values'
        res_str=res_str_stored
    endif
endif

if(grism eq 'H' and filter eq 'H' and n_elements(res_str_stored) eq 1) then begin
    u=[0,1,4,5,7,9,10,11,12,13,14,15]
    c2kywrap_inv=total(((res_str.kywrap_inv/res_str_stored.kywrap_inv)[u]-1)^2)
    if(c2kywrap_inv gt 2.0 and (not keyword_set(force))) then begin
        message,/inf,'Problem computing distorion map. Using stored values'
        res_str=res_str_stored
    endif
endif

if(grism eq 'H' and filter eq 'H3000' and n_elements(res_str_stored) eq 1) then begin
    u=[0,1,4,5,7,9,10,11,12,13,14,15]
    c2kywrap_inv=total(((res_str.kywrap_inv/res_str_stored.kywrap_inv)[u]-1)^2)
    if(c2kywrap_inv gt 2.0 and (not keyword_set(force))) then begin
        message,/inf,'Problem computing distorion map. Using stored values'
        res_str=res_str_stored
    endif
endif

if(grism eq 'K' and filter eq 'K3000' and n_elements(res_str_stored) eq 1) then begin
    u=[0,1,4,5,7,9,10,11,12,13,14,15]
    c2kywrap_inv=total(((res_str.kywrap_inv/res_str_stored.kywrap_inv)[u]-1)^2)
    if(c2kywrap_inv gt 2.0 and (not keyword_set(force))) then begin
        message,/inf,'Problem computing distorion map. Using stored values'
        res_str=res_str_stored
    endif
endif

if(grism eq 'J' and filter eq 'J' and n_elements(res_str_stored) eq 1) then begin
    u=[0,1,4,5,7,9,10,11,12,13,14,15]
    res_str_stored.mask_dy0=-20.0
    c2kywrap_inv=total(((res_str.kywrap_inv/res_str_stored.kywrap_inv)[u]-1)^2)
    if(c2kywrap_inv gt 2.0 and (not keyword_set(force))) then begin
        message,/inf,'Problem computing distorion map. Using stored values'
        res_str=res_str_stored
    endif
endif

slit_reg=get_slit_region(mask,nx=nx,ny=ny,dist_map=res_str,slit_trace=slit_trace)

slits_wlsol = mmirs_wlsol_iguess(grism,mask.x,period=period)
slits_wstart = slits_wlsol[0,*]
slits_wend = slits_wlsol[0,*]+slits_wlsol[1,*]*Nx
prof_mod=dblarr(Ny)
for i=0,n_slit-1 do begin
    n_min=(slit_trace[i].y_trace[Nx-(xref_min+xref_max)/2]-slit_trace[i].slit_h/2.0) > 0
    n_max=(slit_trace[i].y_trace[Nx-(xref_min+xref_max)/2]+slit_trace[i].slit_h/2.0) < (Ny-1)
    wl_ref_min=mask[i].wstart+(slits_wend[i]-slits_wstart[i])*xref_min/2048.0
    wl_ref_max=mask[i].wstart+(slits_wend[i]-slits_wstart[i])*xref_max/2048.0
    if(wl_ref_max ge lin_params.wl_min and wl_ref_min le lin_params.wl_max) then prof_mod[n_min:n_max]=1.0
endfor
psf_krnl=psf_gaussian(ndim=1,fwhm=[2.0],npix=9,/norm)
prof_mod=convol(prof_mod,psf_krnl)

nsegments=4
overlap=0.3
z=5.0
dy=-20.0+findgen(40*z+1L)/z

dy_shift=dblarr(nsegments)
c_shift=dy_shift
y_grid_shift=dblarr(nsegments)

for i=0,nsegments-1 do begin
    ymin = (i eq 0)? skippix : (double(i)-overlap)*Ny/nsegments
    ymax = (i eq nsegments-1)? Ny-skippix-1 : (double(i)+1.0+overlap)*Ny/nsegments
    if(ymin lt 0) then ymin=0
    if(ymax ge Ny) then ymax=Ny-1
    ymin=fix(ymin)
    ymax=fix(ymax)
    y_grid_shift[i]=(ymax+ymin)/2.0
    vec1=rebin(prof_ref[ymin:ymax],(ymax-ymin+1L)*z)
    vec2=rebin(prof_mod[ymin:ymax],(ymax-ymin+1L)*z)
    cf=c_correlate(vec1,vec2,dy*z)
    nn=max(cf,cnn)
    if(cnn eq 0 or cnn eq (n_elements(dy)-1L) or nn lt 0.03) then begin
        dy_shift[i]=!values.f_nan
        c_shift[i]=!values.f_nan
    endif else begin
        dy_shift[i]=dy[cnn]
        c_shift[i]=nn
    endelse
endfor

;coeff_shift=poly_fit(y_grid_shift,dy_shift,1)
;res_str.mask_dy0=res_str.mask_dy0+coeff_shift[0]
;res_str.mask_y_scl=(res_str.mask_y_scl)*(1d -coeff_shift[1])

good_shift=where(finite(dy_shift) eq 1,cgood_shift)
if(cgood_shift eq 0) then dy_shift[*]=0.0
if(cgood_shift gt 0 and cgood_shift le 3) then dy_shift[*]=median(dy_shift[good_shift])

coeff_shift=poly_fit(Ny-y_grid_shift[where(finite(good_shift) eq 1)],dy_shift[where(finite(good_shift) eq 1)],1)
res_str.mask_dy0=res_str.mask_dy0-coeff_shift[0]
res_str.mask_y_scl=(res_str.mask_y_scl)*(1d +coeff_shift[1])

;read,aaa
;stop
writefits,wdir+'dist_map.fits',0
mwrfits,res_str,wdir+'dist_map.fits'
end
