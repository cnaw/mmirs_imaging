function remove_arc_lines,lines_arc,skip_lines
   w=0.1 ;;; 0.1A
   lines_wl=transpose(lines_arc[0,*])
   flag=byte(lines_wl*0)
   if(n_params() eq 1) then skip_lines=[-1.0]
   for i=0,n_elements(skip_lines)-1 do begin 
       b_l=where(abs(lines_wl-skip_lines[i]) lt w, cb_l)
       if(cb_l gt 0) then flag[b_l]=1
   endfor
   return,where(flag eq 0)
end

function ident_arc,wdir,arc_image,lambda0,linetab,pos_lines,lines,fwhm,$
    titl,head_image=head_image,skip_lines=skip_lines,wlmin=wlmin,wlmax=wlmax,$
    plot=plot,logfile=logfile,intens=intens,lines_arc=lines_arc,crossfit=crossfit,smy=smy,$
    mean_curve=mean_curve,pssuff=pssuff,obsflux_lines=obsflux_lines,oh=oh

;arc line spectrum identification
;arc_image		arc line image
;lambda0		initial guess for the wavelength scale
;linetab		arc line list filename
;pos			2D-array: position of lines
;lines			1D-array: wavelengths of lines
;titl			title of output plot

lambda=lambda0
status=0
if(n_elements(pssuff) ne 1) then pssuff=''

if(n_elements(smy) ne 1) then smy=0 ;;; smooth along Y direction, used only if /crossfit is set
arc_image_cur=arc_image
instrument=def_inst(logfile)
slit_id=def_slit(logfile)
!P.background=2^24-1 & !P.color=0
sw_g=0
if keyword_set(plot) and plot eq 'plot' then sw_g=1
a=size(arc_image_cur)
;;
;; cnaw 2017-07-14
;; for some reason arc_image is a 1-D array
n_dimensions = a[0]
Nx=a[1]
Ny=a[2]
if(n_dimensions eq 1) then Ny = 1

print,'size(arc_image_cur) pre transpose ', size(arc_image_cur)
if(smy gt 0) then begin
    for i=0,Nx-1 do arc_image_cur[i,*]=transpose(median(transpose(arc_image_cur[i,*]),smy))
endif
print,'size(arc_image_cur) after transpose ', size(arc_image_cur)
w=1
tresh=3
;; cnaw 2017-07-14
;; idl8 thinks this is a lambda function:
;;d_lambda=(lambda(n_elements(lambda)-1)-lambda(0))/(n_elements(lambda)-1)
;;
d_lambda=(lambda[n_elements(lambda)-1]-lambda[0])/(n_elements(lambda)-1)
x=findgen(Nx)
xflag=bytarr(Nx)

;background subtraction
Ndeg=5
bad_arc=where(finite(arc_image_cur) eq 0,cbad_arc) 

print,'ident_arc: Nx, Ny, cbad_arc, bad_arc ', nx, ny, cbad_arc, bad_arc

if (cbad_arc gt 0) then arc_image_cur[bad_arc]=0
for y=0,Ny-1 do begin
   print,'ident_arc: for y=0,ny-1 loop: y = ', y,' Nx Ny ', nx, ny, size(arc_image_cur)
    f=goodpoly(double(x),double(arc_image_cur[*,y]),Ndeg,1,bgr)
    arc_image_cur[*,y]=arc_image_cur[*,y]-bgr
;;    vector=arc_image_cur
;;    find_peak,x,vector,0,ipix,xpk,ypk,bkpk,ipk
;;    xflag[xpk]=2
;;    xflag=byte(convol(double(xflag),dblarr(5)+1.0))
;;    gflag=where(xflag eq 0)
;;    f=goodpoly(double(x[gflag]),double(arc_image_cur[gflag,y]),Ndeg,1)
;;    arc_image_cur[*,y]=arc_image_cur[*,y]-poly(x,f)
endfor

robomean,arc_image_cur,3,0.5,img_level,img_rms
neg_img=where(arc_image_cur lt -5*img_rms,cneg_img)
if(cneg_img gt 0) then arc_image_cur[neg_img]=0.0

;; cnaw 2017-07-14
;; creating a reference spectrum
;; if it is 2D
if(Ny gt 1) then begin
   arc_obs=median(arc_image_cur[*,Ny/2-2:Ny/2+2],dim=2)
endif else begin
   arc_obs=arc_image_cur
endelse

;; determining the curvature of arc lines as a function of wavelength by making
;; a piece-wise cross-correlation with the reference spectrum

M = 100
x_cross=findgen(2*M+1)-M
mean_curve=fltarr(Ny)
n_seg=5
overlap=0.3
skippix=0.0
curve_arr=fltarr(n_seg,Ny)
curve_img=fltarr(Nx,Ny)

for L=0,Ny-1 do begin
        x_curve=fltarr(n_seg)
        for i=0,n_seg-1 do begin
            nmin=(i eq 0)? skippix : (double(i)-overlap)*Nx/n_seg
            nmax=(i eq n_seg-1)? Nx-skippix-1 : (double(i)+1.0+overlap)*Nx/n_seg
            if(nmin lt 0) then nmin=0
            if(nmax ge Nx) then nmax=Nx-1L
            x_curve[i]=(nmax+nmin)/2.0
            cross_cur=reverse(c_correlate($
                           abs(arc_image_cur[nmin:nmax,L])*COSIN_APOD((nmax-nmin+1),10),$
                           abs(arc_obs[nmin:nmax])*COSIN_APOD((nmax-nmin+1),10),x_cross))
            if(max(cross_cur,/nan) ge 0.3) then begin
                gau=gaussfit(x_cross,cross_cur,G,nterms=3)
                curve_arr[i,L]=(abs(G[1]) lt max(x_cross)*0.8)? G[1] : !values.f_nan
            endif else curve_arr[i,L]=!values.f_nan
        endfor
endfor

ccs=sfit_2deg(curve_arr,3,3,xgrid=x_curve,ygrid=findgen(Ny),kx=kx)
good_pnts=where(finite(curve_arr) eq 1)
dcurve=robust_sigma((curve_arr-ccs)[good_pnts])
bcurve=where(abs(curve_arr-ccs) gt 3*dcurve,cbcurve)
if(cbcurve gt 0) then begin
    curve_arr_tmp=curve_arr
    curve_arr_tmp[bcurve]=!values.f_nan
    ccs=sfit_2deg(curve_arr_tmp,3,3,xgrid=x_curve,ygrid=findgen(Ny),kx=kx)

    dcurve=robust_sigma((curve_arr-ccs)[good_pnts])
    bcurve=where(abs(curve_arr-ccs) gt 3*dcurve,cbcurve)
    if(cbcurve gt 0) then begin
        curve_arr_tmp=curve_arr
        curve_arr_tmp[bcurve]=!values.f_nan
        ccs=sfit_2deg(curve_arr_tmp,3,3,xgrid=x_curve,ygrid=findgen(Ny),kx=kx)
    endif
endif
curve_img=poly2d(findgen(Nx),findgen(Ny),kx,deg1=3,deg2=3)
mean_curve=curve_img[Nx/2,*]

;endif
y=findgen(Ny)
max_curve=abs(min(mean_curve,/nan)-max(mean_curve,/nan))

if sw_g eq 1 then begin
    set_plot,'x'
    window,2,xsize=1000,ysize=550
    !P.background=2^24-1
    !P.color=0
endif
if (sw_g eq 0) then begin
    set_plot,'ps'
    device,file=wdir+'arc'+pssuff+'.ps',/landscape
endif

;;loadct,0 ;;; greyscale colour table
plot,[0,Nx-1],[0,Ny-1],/nodata,xst=1,yst=1,position=[0,0.6,1,1],xcharsize=1E-10,/norm
Npix_Y=(sw_g eq 1)? 220 : Ny
tv,255-bytscl(congrid(arc_image_cur,1000,Npix_Y),img_level-img_rms,img_level+20*img_rms),0,0.6,xsize=1,ysize=0.4,/norm

;;; finding emission lines, starting from the slit centre, going up then down the slit
pos_lines_tmp=dblarr(Ny,5000)
pos_lines_all=dblarr(Ny,5000)
obsflux_lines_all=dblarr(Ny,5000)
pos_number=intarr(Ny)
kzz=1.0 ;;; 1.0
oversample=5.0
kzz=kzz*oversample
for k=Ny/2,Ny-1 do begin
    vector=(arc_image_cur[*,k] > 0)
    find_peak,x,vector,0,ipix,xpk,ypk,bkpk,ipk 
    s_int=reverse(sort(ypk))
    ;; centroid computation
    ;cntrd_1d,vector,xpk,xpk_cnt,fwhm*kzz,/silent,oversample=1
    gcntrd_1d,congrid(vector,Nx*oversample,cub=-0.5),xpk*oversample,xpk_cnt,fwhm*kzz,/silent
    ;cntrd_1d_simple,vector,xpk,xpk_cnt,fwhm*kzz*2.0,/silent,oversample=1
    ypk=ypk[s_int]
    xpk=xpk[s_int]
    xpk_cnt=xpk_cnt[s_int]/oversample
    pos_lines_tmp[k,0:ipk-1]=xpk_cnt-curve_img[(0 > xpk_cnt < (Nx-1)),k]
    pos_lines_all[k,0:ipk-1]=xpk_cnt
    obsflux_lines_all[k,0:ipk-1]=ypk
    if(k eq Ny/2) then intens=ypk
    pos_number[k]=ipk
endfor
for k=0,Ny/2-1 do begin
    vector=(arc_image_cur[*,k] > 0)
    find_peak,x,vector,0,ipix,xpk,ypk,bkpk,ipk
    s_int=reverse(sort(ypk))
    ;cntrd_1d,vector,xpk,xpk_cnt,fwhm*kzz,/silent,oversample=1
    gcntrd_1d,congrid(vector,Nx*oversample,cub=-0.5),xpk*oversample,xpk_cnt,fwhm*kzz,/silent
    ;cntrd_1d_simple,vector,xpk,xpk_cnt,fwhm*kzz*2.0,/silent,oversample=1
    ypk=ypk[s_int]
    xpk=xpk[s_int]
    xpk_cnt=xpk_cnt[s_int]/oversample
    pos_lines_tmp[k,0:ipk-1]=xpk_cnt-curve_img[(0 > xpk_cnt < (Nx-1)),k]
    pos_lines_all[k,0:ipk-1]=xpk_cnt
    obsflux_lines_all[k,0:ipk-1]=ypk
    pos_number[k]=ipk
endfor

;; sorting line positions
print,'Max(pos_number)=',max(pos_number)
w=fwhm
good_pos=fltarr(400)
good_pos_all=fltarr(Ny,400)
j=0
d_Ny=0.33*Ny ;;; maximal number of missing pixels along the slit (33%)
;;d_Ny=0.25*Ny ;;; max.num of missing pixels along the slit

;; finding the nearest peak to the position of every emission line in the reference
;; spectrum taking into account the curvature information
for k=0,max(pos_number)do begin
    pos_line_cur=pos_lines_tmp(Ny/2,k)
    r=where(abs(pos_lines_tmp-pos_line_cur) lt w, index)
    if ((index gt Ny-d_Ny) and (pos_line_cur gt 0)) then begin
        good_pos[j]=pos_line_cur
        good_pos_all[*,j]=pos_lines_all(*,k)
        ;print,'adding k,n_r',k,index,pos_line_cur
        j=j+1
    endif
endfor
r=where(good_pos gt 0,N_line)
good_pos=good_pos[r]
good_pos_all=good_pos_all[*,r]
print,'Number of accepted lines: ',N_line

CASE instrument OF
    'MMIRS': max_N_line = 100
ENDCASE

if(N_line gt max_N_line) then begin
    N_line=max_N_line
    good_pos=good_pos[0:max_N_line-1]
    good_pos_all=good_pos_all[*,0:max_N_line-1]
endif

;; plotting detected lines
for k=0,N_line-1 do begin
        oplot,curve_img[good_pos[k],*]+good_pos[k]-w,findgen(Ny) ,linestyle=1
        oplot,curve_img[good_pos[k],*]+good_pos[k]+w,findgen(Ny) ,linestyle=1
endfor


;; creating the arc model from the initial approximation of the wavelength solution
arc_mod_nn=arc_model(lambda,fwhm,linetab=linetab)
arc_mod=arc_mod_nn/max(arc_mod_nn)

;; offset determination between arc_obs & arc_mod using cross-correlation
C=5.0
M=100.0
x_cross=findgen(2*M*C+1)/C-M
n_nm=n_elements(arc_mod)
n_ns=n_elements(arc_obs)
cross=reverse(c_correlate(congrid(arc_mod,n_nm*C,cubic=-0.5)*COSIN_APOD(Nx*C,20),$
                          congrid(arc_obs,n_ns*C,cubic=-0.5)*COSIN_APOD(Nx*C,20),x_cross))
if(keyword_set(oh)) then begin
    gau=gaussfit(x_cross,cross,G)
    message,/inf,'Wavelength solution guess needs to be adjusted by '+string(G[1],format='(f6.1)')+' pixels'
    if(abs(G[1]) lt 5 and slit_id ne 'mos') then begin
        lambda=lambda+G[1]*d_lambda
    endif else message,/inf,'Ignoring wavelength solution initial guess adjustment'
endif
m_corr=max(cross,m_idx)
print,'dlam=',d_lambda*x_cross[m_idx],m_idx,n_elements(cross)
lambda=lambda+d_lambda*x_cross[m_idx]/C

;; reading the line list and removing those in the "SKIPLINE" keyword
CASE instrument OF
    'MMIRS': BEGIN
                lines_arc=read_asc(linetab)
                if(keyword_set(oh)) then lines_arc=filter_linelist(lines_arc,fwhm*d_lambda) ;;;; compute effective wavelength for OH lines
                if(n_elements(wlmin) eq 1) then lines_arc=lines_arc[*,where(lines_arc[0,*] ge wlmin)]
                if(n_elements(wlmax) eq 1) then lines_arc=lines_arc[*,where(lines_arc[0,*] le wlmax)]

                if(n_elements(skip_lines) eq 0) then begin
                    skip_lines=[-1.0]
                endif
                lines_idx=remove_arc_lines(lines_arc,skip_lines)

                lines=float(transpose(lines_arc[0,lines_idx]))
                index_line=transpose(lines_arc[1,lines_idx])
                flux_line=transpose(lines_arc[2,lines_idx])
                lines_arc=lines_arc[*,lines_idx]
                s_flux=reverse(sort(flux_line))
                lines=lines[s_flux]
                index_line=index_line[s_flux]
                flux_line=flux_line[s_flux]
            ENDCASE
ENDCASE
print,' max_curve',max_curve

;; extracting arc lines from the line list in the spectral range
left=lambda[0]-4*fwhm*d_lambda
right=lambda[Nx-1]-max_curve+(4*fwhm)*d_lambda
wlsel=where(lines gt left and lines lt right)
lines_all=lines[wlsel]
index_line_all=index_line[wlsel]
flux_line_all=flux_line[wlsel]
N_line=N_elements(lines_all)

good_pos_orig=good_pos

;; performing three iterations to identify arc lines

for iter=0,2 do begin
    index=3-iter

    indsel=where(index_line_all ge index)
    lines=lines_all[indsel]
    index_line=index_line_all[indsel]
    flux_line=flux_line_all[indsel]

    used_lines=lines*0+1
    good_line=fltarr(N_elements(good_pos_orig))

    for k=0,N_elements(good_pos_orig)-1 do begin
        n_fwhm=2.0 ;;; was 4.0
        if(iter gt 0) then n_fwhm=1.0/iter
        r=where((abs(lines-lambda[good_pos_orig[k]]) lt n_fwhm*fwhm*d_lambda) and $
                 (used_lines eq 1),ind)
        if (ind eq 1) then begin
            good_line[k]=lines[r]
            used_lines[r]=0
        endif
        if (ind gt 1) then begin ;; more than one line found around the position
            print,'>1 lines: lam,linepos=',lambda[good_pos_orig[k]],lines[r]
            if(iter gt 0) then begin
                mm=min(abs((lines-lambda[good_pos_orig[k]])[r]),midx)
                good_line[k]=lines[r[midx]]
                used_lines[r[midx]]=0
            endif else begin
                r_add=where(flux_line[r]/flux_line[r[0]] gt 0.75, cra)
                if(cra eq 1) then begin
                    good_line[k]=lines[r[r_add]]
                    used_lines[r[r_add]]=0
                    print,'>1 lines: chosenA: ',lines[r[r_add]]
                endif else begin 
                    mm=min(abs((lines-lambda[good_pos_orig[k]])[r[r_add]]),midx)
                    good_line[k]=lines[r[r_add[midx]]]
                    used_lines[r[r_add[midx]]]=0
                    print,'>1 lines: chosenB: ',lines[r[r_add[midx]]]
                endelse
            endelse
        endif
    endfor
    r=where(good_line gt 0,N_line)
    good_line=good_line[r]
    good_pos=good_pos_orig[r]
    print,N_elements(good_line),' lines identified after iteration ',iter
    if(n_elements(good_line) lt 5) then begin
        message,/inf,'Warning: wavelength solution failed, too few lines could be identified'
        if (sw_g eq 0) then begin
            device,/close
        endif
        return,2
    endif

    ck = robust_poly_fit(good_pos,good_line,3)
    lambda=poly(findgen(Nx),ck)
endfor

;; plotting
out=ALOG(arc_obs/max(arc_obs)+0.1)
out=out-min(out(Nx/2-Nx/4:Nx/2+Nx/4))
out=out/max(out(Nx/2-Nx/4:Nx/2+Nx/4))
plot,out,xst=1,position=[0,0,1,0.55],/norm,/noerase,yst=1,yrange=[0,1.19],$
	title=titl
for k=0,N_line-1 do begin
    oplot, [1,1]*good_pos[k],[0,1],linestyle=1
    xyouts,good_pos[k],1,string(good_line[k],format='(F8.2)'),/data,charsize=0.75,orientation=90
endfor
if (sw_g eq 0) then begin
    device,/close
    ;set_plot,'X'
endif
;; done with the plot

pos_lines=fltarr(Ny,N_line)
obsflux_lines=fltarr(Ny,N_line)
for y=0,Ny-1 do $
    pos_lines(y,*)=good_pos[*]+curve_img[(0 > good_pos < (Nx-1)),y]

lines=good_line
l_sort=sort(lines)
lines=lines[l_sort]
pos_lines=pos_lines[*,l_sort]
obsflux_lines=obsflux_lines[*,l_sort]

print,'Filling pos_lines...',format='(a,$)'
for y=0,Ny-1 do begin
    for j=0,N_line-1 do begin
        m=min(abs(pos_lines_all[y,*]-pos_lines[y,j]),cm)
        pos_lines[y,j]=(m lt 30)? pos_lines_all[y,cm] : !values.f_nan
        obsflux_lines[y,j]=(m lt 30)? obsflux_lines_all[y,cm] : !values.f_nan
    endfor
endfor
print,'done'

return, status
end
