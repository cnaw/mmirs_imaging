function crea_disper_mmirs_slit,logfile,slit,fwhm,N_deg,xpos=xpos,ypos=ypos,YFIT=yfit,oh=oh,$
    PLOT=plot,NdegY=NdegY,smy=smy,pssuff=pssuff,ymin_out=ymin_out,ymax_out=ymax_out
;+
; NAME:
;	CREA_DISPER
; PURPOSE:
;      Identification of arc/sky lines and wavelength solution computation for one slit
; DESCRIPTION:
;	TBW
;
; CALLING SEQUENCE:
;	status=CREA_DISPER_MMIRS_SLIT(logfile,slit,fwhm,N_deg,PLOT=plot,$
;                                     yfit=yfit,NdegY=NdegY)
;
; CATEGORY:
;	CfA MMIRS pipeline
;
; INPUTS:
;	LOG = String scalar of file name  FITS-header from the LOG observation
;       SLIT = number of the slit to work on
;	FWHM = spectral line expected fwhm in pixels
;	N_deg = degree of the polynomial along dispersion
;
; OUTPUTS:
;	 2-dimensional ( N_deg+1) x N_spectra float array coefficient dispersion curve,
;        saved to the disk in  working directory with  standard name 'disper.fits' and
;	 print in file 'disper.txt'.
;	 Last value in every string is rms approximation i px
;
; OPTIONAL OUTPUT:
;	no
;
; OPTIONAL INPUT KEYWORDS:
;	PLOT - if present 2-d errors of approximation plotted to display,
;	else save plot to POSTSCRIPT file in working directory with name 'err_line.ps'
;
; RESTRICTIONS:
;	no
; NOTES:
;	no
;
; PROCEDURES USED:
;	Functions :  TBW
;	Procedures:  TBW
;
; MODIFICATION HISTORY:
;       Written by Igor Chilingarian, CfA, Nov 2011
;-
;

status=0

if not(keyword_set(N_deg)) then N_deg=3
if(n_elements(NdegY) ne 1) then NdegY=3
message,'Creating wavelength solution...',/cont
wdir=def_wdir(logfile)
instrument=def_inst(logfile)
grism=def_grism(logfile,filter=filter)
sw_g=0
sw_f=0
if keyword_set(plot) then sw_g=1
if keyword_set(yfit) and yfit eq 'Yfit' then sw_f=1

if(n_elements(pssuff) ne 1) then pssuff=''

;;;; x position of the slit in the mask
if(n_elements(xpos) ne 1) then xpos=0.000
if(n_elements(ypos) ne 1) then ypos=0.000


base_type = (keyword_set(oh))? 'obj' : 'arc'
;; flatfielded and dark-subtracted arc
arc_full_t=readfits(wdir+base_type+'_ff.fits',h_full)*0.0

;; Find semester data were taken
period=get_mmirs_period(wdir+base_type+'_ff.fits',ext=0)

;; Read distortion map
dist_map=mrdfits(wdir+'dist_map.fits',1,/silent)

arc_image_name=wdir+base_type+'_slits.fits'
print,'crea_disper_mmirs_slit : arc_image_name is ', arc_image_name

arc_image=mrdfits(arc_image_name,slit+1,h,/silent)

y_off=sxpar(h,'YOFFSET')
ys_cur=sxpar(h,'NAXIS2')
print,'crea_disper_mmirs_slit: SIZE(arc_image)',size(arc_image)
;;; cosmic ray and hot pixel cleaning -- absolutely required when using night sky to compute the wavelength solution
gain=2.1
rdnoise=12.0
noise_arc = sqrt(djs_median(arc_image,width=3,boundary='reflect')*sxpar(h_full,'EXPTIME')*gain+2*rdnoise^2)/sxpar(h_full,'EXPTIME')/gain
la_cosmic_array,arc_image,noise=noise_arc,gain=gain,readn=rdnoise,maskarr=arc_mask,sigclip=(grism eq 'HK' or grism eq 'H')? 18.0 : 10.0
bpix_arc=where(arc_mask ne 0,cbpix_arc)
if(cbpix_arc gt 0) then begin
    arc_image_med=arc_image*0.0
    arc_med_prof=median(arc_image,dim=2)
    for i=0,n_elements(arc_image[0,*])-1 do arc_image_med[*,i]=arc_med_prof
    arc_image[bpix_arc]=arc_image_med[bpix_arc]
endif
;;; end of cosmic ray cleaning
;;print,'crea_disper_mmirs_slit: after LA_cosmic: y_off, ys_cur', y_off, ys_cur
;;print,'crea_disper_mmirs_slit: size(arc_full) ',size(arc_full)
;;print,'crea_disper_mmirs_slit: size(arc_image)',size(arc_image)

extract_slit_image,arc_full_t,arc_image,y_off,/fill
;;print,'crea_disper_mmirs_slit: post extract_slit_image: size(arc_image)',size(arc_image)
;;mwrfits,arc_image,'/data1/mmirs/reduced/2010.1030/slit_pre_rot.fits',/create
;;
;; rectify slit using distortion map 
;;
arc_r=rotate(poly_2d(rotate(arc_full_t,5),dist_map.kxwrap,dist_map.kywrap,2,cubic=-0.5),5)
extract_slit_image,arc_r,arc_image,y_off,ys_cur
;;print,'crea_disper_mmirs_slit: post extract_slit_image: size(arc_image)',size(arc_image)
;;mwrfits,arc_image,'/data1/mmirs/reduced/2010.1030/slit_post_rot.fits',/create

;;print,'crea_disper_mmirs_slit: post extract_slit_image: size(arc_r)    ',size(arc_r)
;;print,'crea_disper_mmirs_slit: post extract_slit_image: size(arc_image)',size(arc_image)

if(keyword_set(oh)) then begin
    realarc_full_t=arc_full_t
    realarc_image=mrdfits(wdir+'arc_slits.fits',slit+1,/silent)
    extract_slit_image,realarc_full_t,realarc_image,y_off,/fill
    arc_r=rotate(poly_2d(rotate(realarc_full_t,5),dist_map.kxwrap,dist_map.kywrap,2,cubic=-0.5),5)
    extract_slit_image,arc_r,realarc_image,y_off,ys_cur
endif

mwrfits,arc_image,'/data1/mmirs/reduced/2010.1030/slit_post_la_cosmic.fits',/create


arc_r=0

;; scruntch spectrum along the wavelength direction
slit_im_prof=total(arc_image,1,/nan)
plot, slit_im_prof
slit_id = def_slit(logfile)
print,'crea_disper_mmirs_slit:slit_id: ',slit_id
if(slit_id ne 'mos') then slit_im_prof=median(slit_im_prof,25)

;; calculate the median of scruntched spectrum for values > 0
med_slit_im_prof=median(slit_im_prof > 0)

if(med_slit_im_prof eq 0.0 and keyword_set(oh)) then begin
    slit_im_prof=total(realarc_image,1,/nan)
    if(slit_id ne 'mos') then slit_im_prof=median(slit_im_prof,25)
    med_slit_im_prof=median(slit_im_prof > 0)
endif

;; find the range of Y that has good values in profile
;; this breaks down if most pixels have negative values
good_slit=where(slit_im_prof gt 0.3*med_slit_im_prof,cgood_slit)
ymin=good_slit[0]
ymax=good_slit[cgood_slit-1]

print,'crea_disper_mmirs_slit:med_slit_im_prof', med_slit_im_prof
print,'crea_disper_mmirs_slit: cgood_slit,ymin, ymax ',cgood_slit, ymin, ymax

ymin_out=ymin+y_off
ymax_out=ymax+y_off

arc_image=arc_image[*,ymin:ymax]
sxaddpar,h,'NAXIS2',n_elements(arc_image[0,*])

Nx=sxpar(h,'NAXIS1')
Ny=sxpar(h,'NAXIS2')
Nyorig=Ny


print,'crea_disper_mmirs_slit: SIZE(arc_image), ymin,ymax ',size(arc_image)

writesuff=''

barc=where(finite(arc_image) ne 1, bcnt)
if(bcnt gt 0) then arc_image[barc]=0

if(keyword_set(oh)) then begin
    realarc_image=realarc_image[*,ymin:ymax]
    brealarc=where(finite(realarc_image) ne 1, brealcnt)
    if(brealcnt gt 0) then realarc_image[brealarc]=0
endif



slitthr=0
bad_slitthr=bytarr(Ny)
good_slitthr_idx=findgen(Ny)

pixmm=0.033198d

l_par = linearisation_params(grism,filter)
wlmin=l_par.wl_min
wlmax=l_par.wl_max

order = (grism eq 'HK' and filter eq 'zJ')? 2 : 1
dis_coeff=mmirs_wlsol_iguess(grism,xpos,ypos=ypos,order=order,period=period)

print,'crea_disper_mmirs_slit: INITIAL GUESS:',dis_coeff
lambda=poly(findgen(Nx),dis_coeff)
red_end=where(lambda gt wlmax,crwl)
if(crwl gt 0) then begin ;;;; cutting red end
    Nx=red_end[0]
    arc_image=arc_image[0:Nx-1,*]
    lambda=lambda[0:Nx-1]
    if(keyword_set(oh)) then realarc_image=realarc_image[0:Nx-1,*]
endif

lambda_c=lambda[Nx/2]
d_lambda=lambda_c-lambda[Nx/2-1]

name_table=(keyword_set(oh))? '~/idl/mmirs/pipeline/calib_MMIRS/linelists/linesOH_R2k.tab' : '~/idl/mmirs/pipeline/calib_MMIRS/linelists/linesArGNIRS.tab'
date=strcompress(sxpar(h,'date-obs'))
A=sxpar(h,'HISTORY')
A = (n_elements(A) ge 3)? A[2] : ' '
CASE instrument OF
    'MMIRS': sptype=(keyword_set(oh))? 'OH' : 'Ar'
ENDCASE

titl='Spectrum '+sptype+' '+date+' grism '+grism+' Slit#'+string(slit+1,format='(i2.2)')

plot_flag = (sw_g eq 1)? 'plot' : 'no'
tablog=readlog(logfile)
skiplkw='SKIPLINE'
val=sxpar(tablog,skiplkw,count=cntval)
if(cntval eq 1) then skip_lines=(strcompress(val,/remove_all) ne '')? double(strsplit(val,',',/extract)) : [-1.0]

status_id=ident_arc(wdir,arc_image,lambda,name_table,pos_lines,lines,fwhm,titl,$
    PLOT=plot_flag,logfile=logfile,intens=intens,wlmin=wlmin,wlmax=wlmax,$
    lines_arc=lines_arc,crossfit=crossfit,smy=smy,skip_lines=skip_lines,$
    head_image=h,mean_curve=mean_curve,pssuff=pssuff,obsflux_lines=obsflux_lines,oh=oh)
if(status_id ne 0 and keyword_set(oh)) then begin
    message,/inf,'Wavelength solution failed in Slit #'+string(slit,format='(i2.2)')+' using OH lines. Trying to use arcs'
    status_id=ident_arc(wdir,realarc_image,lambda,'~/idl/mmirs/pipeline/calib_MMIRS/linelists/linesArGNIRS.tab',pos_lines,lines,fwhm,titl,$
        PLOT=plot_flag,logfile=logfile,intens=intens,wlmin=wlmin,wlmax=wlmax,$
        lines_arc=lines_arc,crossfit=crossfit,smy=smy,skip_lines=skip_lines,$
        head_image=h,mean_curve=mean_curve,pssuff=pssuff,obsflux_lines=obsflux_lines)
    if(status_id ne 0) then begin
        message,/inf,'Wavelength solution failed in Slit #'+string(slit,format='(i2.2)')+' using arc lines. Returning NaN'
        return,2
    endif else arc_image=realarc_image
endif

N_line=N_elements(lines)

w=4
ycrd=dindgen(Ny)
crossfit=(sw_f eq 1)? 1 : 0

itcur=0
itmax=3
cbad_pos=1
inpbadflag=intarr(n_line)

if(crossfit eq 1) then begin
    n_sig=3.0
    writefits,wdir+'pos_lines_orig.fits',pos_lines
    writefits,wdir+'obsflux_lines_orig.fits',obsflux_lines
    usable_lines=where(lines_arc[1,*] gt 0)
    badflag=bytarr(n_elements(lines))
    inpbadflag=badflag

    cbad_pos=1
    itcur=0
    itmax=10
    while(itcur lt itmax and cbad_pos gt 0) do begin
        itcur++
        print,''
        print,'Iteration ',itcur
        gf=where(badflag eq 0, cgf)
        if(cgf eq 0) then message,'Dispersion relation cannot be built, all arc lines are marked bad'

        kcfit=robust_poly_fit(pos_lines[Ny/2,gf],lines[gf],(N_deg < 5), ykcfit)
        dev_arr=lines[gf]-ykcfit
        d_stdev=robust_sigma(dev_arr)
        print,'Stdev(wave.sol.)=',d_stdev
        bad_pos = where(abs(lines[gf]-ykcfit) gt n_sig*d_stdev, cbad_pos)
        cbad_pos0=cbad_pos
        for j=0,cbad_pos0-1 do begin
            min_dlam=min(abs(lines_arc[0,usable_lines]-(lines[bad_pos[j]]+dev_arr[bad_pos[j]])),mdlidx)
            print,'n_line,min_dlam,newlam=',bad_pos[j],min_dlam,lines_arc[0,usable_lines[mdlidx]],(lines[bad_pos[j]]+dev_arr[bad_pos[j]])
            if(min_dlam lt 15.0*2) then begin
                lines[bad_pos[j]]=lines_arc[0,usable_lines[mdlidx]]
                badflag[bad_pos[j]]=0
                cbad_pos=cbad_pos-1
;;;;                gf=[gf,bad_pos[j]]
                print,'Identified as '+string(lines_arc[0,usable_lines[mdlidx]])
            endif else begin
                badflag[bad_pos[j]]=1
                print,'Left unidentified'
            endelse
        endfor
    endwhile
    gf=where(badflag eq 0, cgf)
    lines=lines[gf]
    pos_lines=pos_lines[*,gf]
    obsflux_lines=obsflux_lines[*,gf]
    n_line=n_elements(gf)
    badflag=bytarr(n_line)
    bflag_pos_lines=bytarr(Ny,n_line)

    yvec=findgen(Ny)
    for i=0,n_line-1 do begin
        gpos_cur=where(finite(pos_lines[*,i]) eq 1,cgpos_cur,compl=badlc,ncompl=cbadlc)
        if(cgpos_cur ge 5) then begin
           print, 'crea_disper_mmirs_slit.pro : NdegY ', ndegy
            kl=robust_poly_fit(yvec[gpos_cur],pos_lines[gpos_cur,i],NdegY)
            fit_pos_line_cur=poly(yvec,kl)
            d_line_cur=pos_lines[*,i]-fit_pos_line_cur
            badlc=where(abs(d_line_cur) gt 1.5,cbadlc) ;;; reject measurements offset by >1.5pix
        endif
        if(cbadlc gt 0) then bflag_pos_lines[badlc,i]=1
    endfor

    itcur=0
    n_pos_lines=n_elements(pos_lines)
    data_arr=dblarr(3,n_pos_lines)
    
    flag_arr=reform(bflag_pos_lines,n_pos_lines)
    bfl=where(flag_arr ne 0, cbfl)
    data_arr[0,*]=reform(pos_lines,1,n_pos_lines)
    err_data_arr=1.0/((reform(obsflux_lines,1,n_pos_lines) > 0.01)^0.25)
    b_err=where((finite(err_data_arr) ne 1), cb_err)
    if(cb_err gt 0) then err_data_arr[b_err]=max(err_data_arr,/nan)

    print,'Rejected ',cbfl,' measurements due to bad positions'
    if(cbfl gt 0) then data_arr[0,bfl]=!values.f_nan
    data_arr[1,*]=reform(dindgen(Ny) # replicate(1.0,n_line),1,n_pos_lines)
    data_arr[2,*]=reform(replicate(1.0,Ny) # lines,1,n_pos_lines)
    linesstd=transpose(data_arr[2,*]*0.0)
    data_line_id=reform(replicate(1l,Ny) # lindgen(n_line),n_pos_lines)
    
    data_arr[0,*]=(data_arr[0,*]-Nx/2d)/(Nx/2d)
    data_arr[1,*]=(data_arr[1,*]-Ny/2d)/(Ny/2d)
    
    gdata=where(finite(total(data_arr,1)) eq 1, cgdata)
    data_arr=data_arr[*,gdata]
    err_data_arr=err_data_arr[*,gdata]

    data_line_id=data_line_id[gdata]
    linesstd=linesstd[gdata]
    print,'Fitting 2D surface using N points; N=',cgdata
    c_out = 1
    while(itcur lt itmax and c_out gt 0) do begin
        itcur++
        bfit=sfit_2deg(data_arr, err=err_data_arr, N_deg, NdegY, kx=kx_f, /irreg)
        resid=data_arr[2,*]-bfit
        for i=0,n_line-1 do begin
            glc=where(data_line_id eq i,cglc)
            if(cglc gt 0) then begin
                linesstd[glc]=robust_sigma(data_arr[2,glc]-bfit[glc])
;                print,'Line, stdev(A):',lines[i],linesstd[glc[0]]
            endif else badflag[i]=1
        endfor
        medlstd=median(linesstd)
        print,'Median(linesstd)=',medlstd

        res_std=robust_sigma(resid)
        outl=where(((abs(resid) gt n_sig*res_std) or $
                    (linesstd gt medlstd*n_sig*2)) and (data_line_id ne 0), c_out, compl=goodfit)
;;        outl=where((abs(resid) gt n_sig*res_std) and data_line_id ne 0, c_out, compl=goodfit)
        print,'Iteration, res_std, c_out=',itcur,res_std,c_out
        gdata=gdata[goodfit]
        data_arr=data_arr[*,goodfit]
        err_data_arr=err_data_arr[*,goodfit]
        data_line_id=data_line_id[goodfit]
        linesstd=linesstd[goodfit]
    endwhile
    resid=data_arr[2,*]-bfit[goodfit]

    mask_posline=fltarr(n_elements(pos_lines))+!values.f_nan
    mask_posline[gdata]=0.0
    pos_lines=pos_lines+mask_posline

    mean_lines_orig=poly2d((reform(pos_lines,n_pos_lines)-Nx/2d)/(Nx/2d),$
                           (reform(dindgen(Ny)#replicate(1.0,n_line),n_pos_lines)-Ny/2d)/(Ny/2d),$
                           kx_f,deg1=N_deg,deg2=NdegY,/irreg)
    mean_lines_orig=reform(mean_lines_orig,Ny,n_pos_lines/Ny)
    mean_lines=mean_lines_orig
    disper_par=dblarr(N_deg+2,Ny)
    yvec=(dindgen(Ny)-Ny/2d)/(Ny/2d)
    pp=-1d + dindgen(6)*0.4d
    pp1=pp*Nx/2d +Nx/2d
    for i=0,Ny-1 do begin
;;        for j=0,N_deg do begin
        for j=0,NdegY do begin
            disper_par[j,i]=poly(yvec[i],kx_f[*,j])
        endfor
        dp=poly_fit(pp1,poly(pp,disper_par[0:N_deg,i]),N_deg)
        disper_par[0:N_deg,i]=dp
        disper_par[N_deg+1,i]=0.00001
    endfor
    pos_lines_orig=pos_lines
    err_pos_lines=pos_lines*0.0+0.001
;    read,aaa
endif else begin
    while(itcur lt itmax and cbad_pos gt 0) do begin
        disper_par=approx_disp_relation(arc_image,lambda,pos_lines,lines,mean_lines,$
                                 N_deg,fwhm,intens=intens,fwhm=fwhm,ndegy=ndegy,smy=smy,$
                                 crossfit=crossfit,ycrd=ycrd,pos_lines_orig=p_ori,$
                                 slitthr=slitthr,badflag=badflag,$
                                 instrument=instrument,inpbadflag=inpbadflag,bad_slitthr=bad_slitthr,err_pos_lines=err_p_ori)
        badlines=where(badflag eq 1, cbadlines, compl=goodlines)
        if(cbadlines gt 0) then print,'Lines rejected: ',lines[badlines],' accepted:',n_elements(goodlines)
        print,'Lines accepted:',n_elements(goodlines)
        if(itcur eq 0) then begin
            pos_lines_orig=p_ori
            err_pos_lines=err_p_ori
            writefits,wdir+'pos_lines_orig'+writesuff+'.fits',pos_lines_orig
        endif

        corr_lines=lines
        ;robomean,disper_par(N_deg+1,*),3,0.5,average_rms
        rms_arr=dblarr(N_line)
        dev_arr=dblarr(N_line)
        mean_lines_orig=pos_lines_orig*0
        for i=0,Ny-1 do mean_lines_orig[i,*]=poly(pos_lines_orig[i,*],disper_par[0:n_deg,i])
    
        for j=0,N_line-1 do begin
            ;;;;robomean,mean_lines_orig(*,j),3,0.5,average
            resistant_mean,mean_lines_orig[good_slitthr_idx,j], 3.0, average
            ;    if ABS(lines(j)-average) gt 0.2*d_lambda then corr_lines(j)=average
            if ABS(lines(j)-average) gt 0.4*d_lambda then corr_lines(j)=average
            def=mean_lines_orig(*,j)-lines(j)
            ;;;;robomean,def,3,0.5,mean_def,rms
            resistant_mean,def[good_slitthr_idx], 3.0, mean_def, rms1, num_rej
            rms=rms1*sqrt(n_elements(good_slitthr_idx)-num_rej)
            ;print,'j,RMS=',lines[j],mean_def,rms
            rms_arr[j]=rms
            dev_arr[j]=mean_def
            ;    read,aaa
        endfor
        max_rms=0.45 ;;; 0.4 was before
        badflag=badflag*0
        inpbadflag=inpbadflag*0
        bad_rms=where(rms_arr gt max_rms, cbad_rms)
        if(cbad_rms gt 0) then begin
            print,'RMS rejected: ',lines[bad_rms],rms_arr[bad_rms]
            inpbadflag[bad_rms]=1
        endif
        bad_pos=where((abs(dev_arr) gt 2.0) and (rms_arr le max_rms), cbad_pos) ;; deviation > 2A
        ;bad_pos=where(abs(dev_arr) gt 2.0, cbad_pos) ;; deviation > 2A
        print,'BAD_POS',bad_pos
        usable_lines=where(lines_arc[1,*] gt 0)
        for j=0,cbad_pos-1 do begin
            min_dlam=min(abs(lines_arc[0,usable_lines]-(lines[bad_pos[j]]+dev_arr[bad_pos[j]])),mdlidx)
            print,'min_dlam,newlam=',min_dlam,lines_arc[0,usable_lines[mdlidx]],(lines[bad_pos[j]]+dev_arr[bad_pos[j]])
            if(min_dlam lt 5.0) then begin
                lines[bad_pos[j]]=lines_arc[0,usable_lines[mdlidx]]
                inpbadflag[bad_pos[j]]=0
                badflag[bad_pos[j]]=0
                print,'Identified'
            endif else begin
                inpbadflag[bad_pos[j]]=1
                badflag[bad_pos[j]]=1
                print,'Left unidentified'
            endelse
        endfor
        itcur=itcur+1
    endwhile
    
    badflag=inpbadflag+badflag
    luniq=uniq(lines,sort(corr_lines))
    lines=lines[luniq]
    corr_lines=corr_lines[luniq]
    pos_lines=pos_lines[*,luniq]
    pos_lines_orig=pos_lines_orig[*,luniq]
    err_pos_lines=err_pos_lines[*,luniq]
    inpbadflag=inpbadflag[luniq]
    badflag=badflag[luniq]
    n_line=n_elements(lines)
    ;print,'INPBADFLAG:',byte(inpbadflag)
    ;print,'BADFLAG:',byte(badflag)
    ;;2-nd approximation dispersion curve
    disper_par=approx_disp_relation(arc_image,lambda,pos_lines,corr_lines,mean_lines,$
                                 N_deg,fwhm,intens=intens,fwhm=fwhm,ndegy=ndegy,smy=smy,$
                                 crossfit=crossfit,ycrd=ycrd,$
                                 slitthr=slitthr,instrument=instrument,bad_slitthr=bad_slitthr,$
                                 inpbadflag=inpbadflag,badflag=badflag)
    badflag=badflag+inpbadflag
    
    disp=total(disper_par(1,*))/Ny
    dif=(total(mean_lines,2)-total(lines))/Ny
endelse


;plot result of approximation
if sw_g eq 1 then window,3,xsize=600,ysize=850
if sw_g eq 0 then begin
    set_plot,'PS'
    device,file=wdir+'err_approx'+pssuff+'.ps',xsize=22,ysize=28,xoffset=0,yoffset=1,/portrait
endif
!P.background=16777215 &   !P.color=0
dy=4 & ymin=-dy & ymax=dy*N_line
plot,[0,Ny-1],[ymin,ymax],xst=1,yst=1,/nodata,$
	ytitle='Residual deviations, px',$
	xtitle='Slit position',$
	position=[0.09,0.05,0.75,0.82],/norm
print,'Lines identified'
print,'Lambda','dLam','RMS','pos','e_pos',format='(a9,a7,a5,2a9)'
openw,u,wdir+'lines_id'+pssuff+'.txt',/get_lun
printf,u,'Lines identified'
printf,u,'Lambda','dLam','RMS','pos','e_pos',format='(a9,a7,a5,2a9)'
mean_lines_orig=pos_lines_orig*0
for i=0,Ny-1 do mean_lines_orig[i,*]=poly(pos_lines_orig[i,*],disper_par[0:n_deg,i])
for j=0,N_line-1 do begin
    def_full=mean_lines_orig[*,j]-lines[j]
    def=mean_lines[*,j]-lines[j]
;;;    aaa='' & read,aaa & if aaa eq 's' then stop
    robomean,def[good_slitthr_idx],3,0.5,mean_def,rms
;;;    resistant_mean,def[good_slitthr_idx],3,mean_def,rms
    oplot,10*def/d_lambda+dy*j,psym=6,symsize=0.2,col=128
    oplot,10*def_full/d_lambda+dy*j,psym=-4,symsize=0.2
    oplot,[0,Ny],[1,1]*dy*j,color=150
    correction=mean_def
    if correction lt 0 then S='  -'+string(abs(correction),format='(F4.2)')
    if correction gt 0 then S='  +'+string(abs(correction),format='(F4.2)')
    if correction eq 0 then S='        '
    flagval=(badflag[j] eq 0)?' ':'*'
    xyouts,Ny+5,dy*j,string(lines[j],format='(F8.2)')+$
        S+string(rms,format='(F7.2)'),charsize=0.75,/data
    print,lines[j],correction,rms,pos_lines[Ny/2,j],err_pos_lines[Ny/2,j],flagval,format='(f9.2,f7.2,f5.2,f9.2,f8.3,1x,a2)'
    printf,u,lines[j],correction,rms,pos_lines[Ny/2,j],err_pos_lines[Ny/2,j],flagval,format='(f9.2,f7.2,f5.2,f9.2,f8.3,1x,a2)'
endfor
close,u
free_lun,u
x_hist=findgen(21)*0.025
plot,x_hist,histogram(disper_par(N_deg+1,*)/d_lambda,min=0,max=.5,binsize=.025),$
	position=[0.09,0.86,0.38,0.98],/noerase,/norm,$
	xst=1,yst=1,psym=10,$
	xtitle='error, px'
xyouts,0.45,0.975,'ACCURACY OF THE 2D WAVELENGTH SOLUTION',/norm
xyouts,0.4,0.95,'date of observation   '+date,/norm
xyouts,0.4,0.93,'Spectrum Ne-Ar-He  file:'+ a,/norm
xyouts,0.4,0.91,'Grating '+strcompress(sxpar(h,'disperse'))+$
	' TILT='+STRING(SXPAR(H,'TILTPOS'),format='(F6.1)'),/norm
cw=0 & disp=total(disper_par(1,*))/Ny
for j=0,N_deg do  cw=cw+disper_par(j,Ny/2)*(Nx/2)^J
xyouts,0.4,0.89,'Central wavelength'+string(cw,format='(F8.2)')+$
	' A Mean dispersion'+string(disp,format='(F6.2)')+' A/px',/norm
;robomean,disper_par(N_deg+1,*),3,0.5,disp_rms
resistant_mean,disper_par(N_deg+1,*),3,disp_rms
xyouts,0.4,0.87,'Number of lines'+string(N_line,format='(I3)')+$
	'   Average error approximation '+string(disp_rms/d_lambda,format='(F4.2)')+' px',/norm
xyouts,0.78,0.82,'!7k!3(A)    !7Dk!3  rms !7k!3',/norm
if sw_g eq 0 then begin
    device,/close
    set_plot,'X'
endif

;saving dispersion relation
writefits,wdir+'disper'+writesuff+'.fits',disper_par
writefits,wdir+'pos_lines'+writesuff+'.fits',pos_lines
a=size(disper_par)
out_format='(F8.2,F8.5,'+string(a(1)-3,format='(I1)')+'E16.6,F6.2)'
openw, UNIT,wdir+'disper'+writesuff+'.txt', /GET_LUN
for k=0,a(2)-1 do begin
    printf, UNIT,disper_par(*,k),format=out_format
endfor
close, UNIT
free_lun,UNIT

message,'dispersion relation mean error ='+string(disp_rms/d_lambda,format='(F5.2)')+' px',/cont

disper_par=0
return,status

end
