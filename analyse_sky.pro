pro analyse_sky,logfile,image_type,raw=raw,plot=plot,$
    bin_s=bin_s,min_v=min_v,max_v=max_v,gain=gain

if(n_params() eq 1) then image_type='obj'

if(n_elements(gain) ne 1) then gain=2.1 ;;;; 5.0
if(n_elements(min_v) ne 1) then min_v=0.0
if(n_elements(max_v) ne 1) then max_v=5.0
if(n_elements(bin_s) ne 1) then bin_s=0.05


log=readlog(logfile)
rdir = sxpar(log,'R_DIR')
wdir = sxpar(log,'W_DIR')

sci_pref = sxpar(log,'SCI')
dithpos = sxpar(log,'DITHPOS')

sci_pref2 = sxpar(log,'SCI2',count=cntval)
diffmode = (cntval eq 1)? 1 : 0

if(keyword_set(raw)) then begin
    rawext=sxpar(log,'RAWEXT',count=cntre)
    if(cntre eq 0) then rawext=''
    f_raw=rdir+get_filename(sci_pref,suffix='.fix.fits'+rawext)
    h_raw=headfits(f_raw)
    raw=mrdfits(f_raw,1,h_raw_ext,/silent)
    sxaddpar,h_raw_ext,'BZERO',0.0
    writefits,wdir+'rawobj_dark.fits',0,h_raw
    mwrfits,float(raw),wdir+'rawobj_dark.fits',h_raw_ext
    
    extract_2d_slits,logfile,'rawobj',/unproc

    pref_name='raw'
endif else pref_name=''

mask=read_mask_mmirs_ms(wdir+'mask_mos.txt')

n_slits = n_elements(mask)

targ_idx=where(mask.type eq 'TARGET',ctarg_idx)

period=get_mmirs_period(wdir+'obj_dark.fits', year=year, exten=1)
detector_id = (period eq '2014B' or year gt 2014)? 'H2RG' : 'H2_56'

if(detector_id eq 'H2RG') then begin
    gain=0.94
endif

good_all_raw=[-1]
good_all_ss_sing=[-1]
good_all_ss_diff=[-1]
for s=0,ctarg_idx-1 do begin
    rawim=mrdfits(wdir+pref_name+'obj_slits.fits',targ_idx[s]+1,/silent)*gain
    ss_sing=mrdfits(wdir+'obj-sky_slits.fits',targ_idx[s]+1,/silent)*gain
    ss_diff=(diffmode eq 1)? $
            mrdfits(wdir+'obj_diff-sky_slits.fits',targ_idx[s]+1,/silent)*gain : $
            ss_sing*0.0
    g_cur=where(finite(rawim+ss_sing+ss_diff) eq 1, cg_cur)

    ;;;; masking the targets should go here
    ;;;;

    if(cg_cur gt 0) then begin
        good_all_raw=[good_all_raw,rawim[g_cur]]
        good_all_ss_sing=[good_all_ss_sing,ss_sing[g_cur]]
        good_all_ss_diff=[good_all_ss_diff,ss_diff[g_cur]]
    endif

endfor

good_all_raw=good_all_raw[1:*]
good_all_ss_sing=good_all_ss_sing[1:*]
good_all_ss_diff=good_all_ss_diff[1:*]

hdr_sing=headfits(wdir+'obj-sky_slits.fits')

exptime=sxpar(hdr_sing,'EXPTIME',count=nval)
if(nval eq 0) then exptime=1.0

rdnoise=sxpar(hdr_sing,'RDNOISE',count=nval)
if(nval eq 0 or rdnoise eq 0.0) then rdnoise=16.0

h_cnt = histogram(alog10(good_all_raw*exptime),$
    min=min_v,max=max_v,binsize=bin_s,reverse_ind=h_cnt_sub)

n_bins=n_elements(h_cnt)

h_ss_sing = dblarr(n_bins)+!values.f_nan
h_ss_diff = h_ss_sing
h_ss_sing_sigcl = h_ss_sing
h_ss_diff_sigcl = h_ss_sing

for i=1,n_bins-1 do begin
    if((h_cnt_sub[i]-h_cnt_sub[i-1]) ge 3) then begin
        h_ss_sing[i-1]=stdev(good_all_ss_sing[h_cnt_sub[h_cnt_sub[i-1]:h_cnt_sub[i]-1L]])
        h_ss_diff[i-1]=stdev(good_all_ss_diff[h_cnt_sub[h_cnt_sub[i-1]:h_cnt_sub[i]-1L]])
        h_ss_sing_sigcl[i-1]=robust_sigma(good_all_ss_sing[h_cnt_sub[h_cnt_sub[i-1]:h_cnt_sub[i]-1L]])
        h_ss_diff_sigcl[i-1]=robust_sigma(good_all_ss_diff[h_cnt_sub[h_cnt_sub[i-1]:h_cnt_sub[i]-1L]])
    endif
endfor

sky_str=replicate({$
                   intens:!values.f_nan,$
                   int_min:!values.f_nan,$
                   int_max:!values.f_nan,$
                   rms_single:!values.f_nan,$
                   rms_diff:!values.f_nan,$
                   rms_single_sigcl:!values.f_nan,$
                   rms_diff_sigcl:!values.f_nan,$
                   numpix:-1L,$
                   t_rms_single:!values.f_nan,$
                   t_rms_diff:!values.f_nan},n_bins)

int_ax=(bin_s*0.5+min_v+findgen(n_bins)*bin_s)

sky_str.intens=10^int_ax*exptime
sky_str.int_min=10^(min_v+findgen(n_bins)*bin_s)*exptime
sky_str.int_max=10^(bin_s+min_v+findgen(n_bins)*bin_s)*exptime
sky_str.rms_single=h_ss_sing*exptime
sky_str.rms_diff=h_ss_diff*exptime
sky_str.rms_single_sigcl=h_ss_sing_sigcl*exptime
sky_str.rms_diff_sigcl=h_ss_diff_sigcl*exptime
sky_str.numpix=h_cnt
sky_str.t_rms_single=sqrt(10^int_ax+rdnoise^2)*exptime
sky_str.t_rms_diff=sqrt(10^int_ax+rdnoise^2)*sqrt(2)*exptime

writefits,wdir+pref_name+image_type+'_sky_sub_stat.fits',0
mwrfits,sky_str,wdir+pref_name+image_type+'_sky_sub_stat.fits'

if(not keyword_set(plot)) then begin
    set_plot,'ps'
    device,/color,filename=wdir+pref_name+image_type+'_sky_sub_stat.ps',xs=18,ys=13
    !x.thick=3
    !y.thick=3
    !p.thick=3
    !p.font=1
endif

loadct,39


ploterror,10^int_ax,h_ss_sing*exptime,10^int_ax/10.0,h_ss_sing*exptime/sqrt(h_cnt),$
        xs=1,psym=-4,/xlog,/ylog,xmar=[6,1],ymar=[3,1],$
        xtitle='Intensity, e-',ytitle='RMS of residuals, e-'
oploterror,10^int_ax,h_ss_diff*exptime,$
        10^int_ax/10.0,h_ss_diff*exptime/sqrt(h_cnt),$
        psym=-4,col=254,errcol=254
oploterror,10^int_ax,h_ss_sing_sigcl*exptime,$
        10^int_ax/10.0,h_ss_sing_sigcl*exptime/sqrt(h_cnt),$
        psym=-7,col=!p.color,errcol=!p.color,linest=1
oploterror,10^int_ax,h_ss_diff_sigcl*exptime,$
        10^int_ax/10.0,h_ss_diff_sigcl*exptime/sqrt(h_cnt),$
        psym=-7,col=254,errcol=254,linest=1
oplot,10^int_ax,sqrt(10^int_ax+rdnoise^2),col=128
oplot,10^int_ax,sqrt(10^int_ax+rdnoise^2)*sqrt(2),col=80
oplot,10^int_ax,h_cnt,col=220,psym=10
;; cnaw 2017-07-19 Since now there is an idl 8.4 procedure with this name,
;; have to change the goddard one to al_legend
al_legend,['Single frame sky subtraction','Difference image sky subtraction',$
        'Theoretical limit for a single frame','Theoretical limit for a difference image',$
        'Histogram of counts'],$
       col=[!p.color,254,128,80,220],psym=[-4,-4,-3,-3,-3],linest=[0,0,0,0,0],box=0

if(not keyword_set(plot)) then begin
    device,/close
    set_plot,disp_family()
    !p.font=-1
    !x.thick=0
    !y.thick=0
    !p.thick=0
endif



end
