;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  MMFIXEN_NONLIN -- IDL implementation of the up-the-ramp fitting routine
;     including the non-linearity correction using pre-computed calibrations
;
;
;
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mmfixen_nonlin,inpfile0,outfile,first=first,linear=linear,keepfirst=keepfirst,$
    verbose=verbose,debug=debug,biasframe=biasframe,badamp=badamp,$
    crosstalk=crosstalk,compress=compress,tmpdir=tmpdir,clean=clean

if(n_elements(tmpdir) ne 1) then tmpdir='/tmp/'
if(n_elements(first) ne 1) then first=2
if(n_elements(badamp) ne 1) then badamp=20 ;;; fixing bad amplifier in the MMIRS H2_56 detector
if(keyword_set(debug)) then verbose=1

tmpfilename=''
if(~file_test(inpfile0)) then begin ;;; no file under the "inpfile" name
    inp_f_arr=file_search(inpfile0+'*')
    if(strlen(inp_f_arr[0]) gt 0) then begin
        inpext=strmid(inp_f_arr[0],strlen(inp_f_arr[0])-3)
        if(inpext eq '.xz') then inpcompress='unxz'
        if(inpext eq '.gz') then inpcompress='gunzip'
        if(inpext eq 'bz2') then inpcompress='bunzip2'
        tmpfilename=tmpdir+'mmirs_pipe_'+string(randomu(seed)*1d+7,format='(i8.8)')+'.fits'
        if(keyword_set(verbose)) then print,'Decompressing '+inp_f_arr[0]+' ...',format='(a,$)'
        spawn,inpcompress+' -c '+inp_f_arr[0]+' >'+tmpfilename
        if(keyword_set(verbose)) then print,'done'
        inpfile=tmpfilename
    endif else begin
        message,'Input file '+inpfile0+' not found. Cannot continue MMFIXEN'
    endelse 
endif else inpfile=inpfile0

h_dummy=headfits(inpfile)

n_img=get_nreads_ramp(inpfile,exptime=exptime,ramp=ramp,header=h_img,gain=gain,rdnoise=rdnoise)
if(rdnoise eq 0.0) then rdnoise=12.0
if(n_img eq 1) then begin
    h_dummy=['SIMPLE  =                    T',$
             'BITPIX  =                   16',$
             'NAXIS   =                    0',$
             'EXTEND  =                    T',$
             'END                           ']
endif

period=get_mmirs_period(inpfile, year=year)
detector_id = (period eq '2014B' or year gt 2014)? 'H2RG' : 'H2_56'
if(detector_id eq 'H2RG') then badamp=-1
if(detector_id eq 'H2RG') then keepfirst=1
if(detector_id eq 'H2RG') then clean=1 ; 1 ;;; perform sigma-clipping for cleaning read-outs

if(n_img le 2) then keepfirst=1 ; force keeping the first readout in the fowler mode

h_img=h_img[where(strlen(strcompress(h_img,/remove_all)) gt 0)]

sxaddpar,h_img,'LONGSTRN','OGIP 1.0',' The OGIP long string convention may be used.',after='XTENSION'

;; read image and correct for reference pixels if applicable
;; for first frame

img=float(mrdfits(inpfile,((n_img eq 1)? 0 : 1),/silent))+32768.0
if(detector_id eq 'H2RG') then img=refpix(img)
nx=long(sxpar(h_img,'NAXIS1'))
ny=long(sxpar(h_img,'NAXIS2'))
sxdelpar,h_img,'BZERO'
sxaddpar,h_img,'BUNIT','counts/s'

if(n_elements(biasframe) eq 1 and detector_id ne 'H2RG') then begin
    if(~file_test(biasframe)) then begin ;;; no file under the "biasframe" name
        bias_f_arr=file_search(biasframe+'*')
        if(strlen(bias_f_arr[0]) gt 0) then begin
            biasext=strmid(bias_f_arr[0],strlen(bias_f_arr[0])-3)
            if(biasext eq '.xz') then biascompress='unxz'
            biasframe_full=bias_f_arr[0]
        endif else begin
            message,'Bias file '+biasframe+' not found. Cannot continue MMFIXEN'
        endelse 
    endif else biasframe_full=biasframe

    bias_n1=3
    bias_n2=4

    if(n_img eq 1) then bias_n1=0
    if(n_img eq 2) then begin ;;; fowler mode
        bias_n1=1
        bias_n2=2
    endif
    if(n_img eq 3) then begin
        bias_n1=3
        bias_n2=1
    endif
    bias = (detector_id eq 'H2RG')? 0.0 : get_bias(biasframe_full,n1=bias_n1,n2=bias_n2,compress=biascompress,h2rg=((detector_id eq 'H2RG')? 1 : 0))
endif else bias=0.0

;; subtract bias and calculate crosstalk correction
;; for first frame

img=img-bias
if(keyword_set(crosstalk)) then $
    img=img-crosstalk_correction_32amp(img)

;; correct for reference pixels, subtract bias, calculate crosstalk
;; for remaining readouts

img_cube=fltarr(nx,ny,n_img)
img_cube[*,*,n_img-(first-1)]=img
if(keyword_set(verbose)) then print,'Reading up-the-ramp cube...',format='(a,$)'
for i=1,n_img-1 do begin
    t=float(mrdfits(inpfile,i+1,/silent))+32768.0-bias
    if(detector_id eq 'H2RG') then t=refpix(t)
    if(keyword_set(crosstalk)) then $
        t=t-crosstalk_correction_32amp(t)
    img_cube[*,*,i-(first-1)]=t
endfor
if(keyword_set(verbose)) then print,'done'
;;;;img_cube=shift(img_cube,0,0,1-first)

;;if(detector_id eq 'H2RG') then refpix_cube,img_cube

if(strlen(tmpfilename) gt 5) then file_delete,tmpfilename,/allow_nonexistent
sxaddpar,h_img,'XTALKCOR',keyword_set(crosstalk),' cross-talk correction status'

;; if linearity correction keyword is not set:

if((not keyword_set(linear)) and (detector_id ne 'H2RG')) then begin
    ;;;; reset effect parameters
    decay_norm=0.0
    decay_time=1.0
    if(ramp eq 1) then begin
        decay_norm = 0.18
        decay_time = 2.7
    endif else begin    ;;;; ramp 5 assumed
        decay_norm = 0.0215*0.0
        decay_time = 2.7
    endelse
    dec_21corr = 1.0 + decay_norm*(exp(-1.0/decay_time)-exp(-2.0/decay_time))

    if(n_img eq 1) then begin
        decay_norm=0.0
        rate_t1=img_cube/dec_21corr
    endif else $
        rate_t1 = (n_img gt 2)? (img_cube[*,*,2]-img_cube[*,*,1])/dec_21corr : (img_cube[*,*,1]-img_cube[*,*,0])/dec_21corr
    if(decay_norm gt 0) then $
        for i=0,n_img-1 do img_cube[*,*,i]=img_cube[*,*,i]-decay_norm*rate_t1*exp(-float(i)/decay_time)
endif

if(not keyword_set(keepfirst)) then begin
    img_cube=img_cube[*,*,1:*]
    n_img=n_img-1
endif

;; seems to be fixing/doing something in the case of bad amplifiers

sat_mask=intarr(nx,ny) ;;; saturation
crh_mask=intarr(nx,ny) ;;; cosmic ray hits
mask_image=bytarr(nx,ny)

sxaddpar,h_img,'SOFTWARE',get_mmirs_pipeline_version()
sxaddpar,h_img,'NREADOUT',n_img,' number of readouts used'

bias_adj=replicate({readout:0,bad_amp:badamp,biasdiff:0.0,used:1},n_img)
bias_adj.readout=findgen(n_img)+(1-keyword_set(keepfirst))

if(badamp ge 0 and n_img ge 3 and ramp ne 7) then begin ;;; excluding ramp=7
    amps_conf=mrdfits('~/idl/mmirs/pipeline/calib_MMIRS/detector/'+detector_id+'/amps32_layout.fits',1,/silent)
    badamp_conf=amps_conf[badamp]
;;; fixing bad amplifier here --- to be written
    if(badamp_conf.direction eq 2) then begin
        n_iter=0
        while(n_iter lt 2) do begin
            fl_xmin=reform(img_cube[badamp_conf.xmin,badamp_conf.ymin:badamp_conf.ymax,*],1024,n_img)
            fl_xmin1=reform(img_cube[badamp_conf.xmin+1,badamp_conf.ymin:badamp_conf.ymax,*],1024,n_img)
            fl_xmax=reform(img_cube[badamp_conf.xmax,badamp_conf.ymin:badamp_conf.ymax,*],1024,n_img)
            fl_xmax1=reform(img_cube[badamp_conf.xmax-1,badamp_conf.ymin:badamp_conf.ymax,*],1024,n_img)
            if(badamp_conf.xmin gt 0 and badamp_conf.xmin ne 1024) then begin
                fl_ref_min=reform(img_cube[badamp_conf.xmin-1,badamp_conf.ymin:badamp_conf.ymax,*],1024,n_img) 
                fl_ref_min1=reform(img_cube[badamp_conf.xmin-2,badamp_conf.ymin:badamp_conf.ymax,*],1024,n_img)
            endif else fl_ref_min=fl_xmin*0.0
            if(badamp_conf.xmax lt 2047 and badamp_conf.xmax ne 1023) then begin
                fl_ref_max=reform(img_cube[badamp_conf.xmax+1,badamp_conf.ymin:badamp_conf.ymax,*],1024,n_img)
                fl_ref_max1=reform(img_cube[badamp_conf.xmax+2,badamp_conf.ymin:badamp_conf.ymax,*],1024,n_img)
            endif else fl_ref_max=fl_xmax*0.0
            dfl_xmin=fl_xmin[*,1:*]-fl_xmin
            dfl_xmin1=fl_xmin1[*,1:*]-fl_xmin1
            dfl_xmax=fl_xmax[*,1:*]-fl_xmax
            dfl_ref_min=fl_ref_min[*,1:*]-fl_ref_min
            dfl_ref_max=fl_ref_max[*,1:*]-fl_ref_max
            
            ;;;; checking difference frames first
            robomean,median((dfl_xmin1-dfl_ref_min)^2,dim=1),3.0,0.5,rmi1,rmi2,rmi3
            robomean,median((dfl_xmax-dfl_ref_max)^2,dim=1),3.0,0.5,rma1,rma2,rma3
            dfl_ref_min_vec=dblarr(n_img-1)
            dfl_ref_max_vec=dblarr(n_img-1)
            for i=0,n_img-2 do begin
                dfl_ref_min_vec[i]=robust_sigma(dfl_ref_min[*,i])
                dfl_ref_max_vec[i]=robust_sigma(dfl_ref_max[*,i])
            endfor
            good_diff = where((median((dfl_xmin1-dfl_ref_min)^2,dim=1) lt 4.0*2.0*dfl_ref_min_vec) or $ ;;; 3 sigma outliers allowed
                              (median((dfl_xmax-dfl_ref_max)^2,dim=1) lt 4.0*2.0*dfl_ref_max_vec),cgood_diff,$
                              compl=bad_diff,ncompl=cbad_diff)
    ;        good_diff = where(abs(median((dfl_xmin-dfl_ref_min)^2,dim=1)) lt 8.0+3.5*sqrt(abs(median(dfl_ref_min))) or $
    ;                          abs(median((dfl_xmax-dfl_ref_max)^2,dim=1)) lt 8.0+3.5*sqrt(abs(median(dfl_ref_max))),cgood_diff,$
    ;                          compl=bad_diff,ncompl=cbad_diff)
    ;;        bad_diff = where(abs(median((dfl_xmin-dfl_ref_min)^2,dim=1)-rmi1) gt 2.0*rmi3 or $
    ;;                         abs(median(dfl_xmax-dfl_ref_max)^2,dim=1)-rma1) gt 2.0*rma3,cbad_diff,$
    ;;                          compl=good_diff,ncompl=cgood_diff)
    
            n_iter++
            bias_diff_arr=((fl_xmin1-(2*fl_ref_min-fl_ref_min1))+$
                              (fl_xmax-(2*fl_ref_max-fl_ref_max1)))/2.0
            bias_diff=median(bias_diff_arr,dim=1)
            ;print,'n_iter=',n_iter,bias_diff,good_diff
            if(cbad_diff eq 0) then begin
                img_cube[badamp_conf.xmin:badamp_conf.xmax,badamp_conf.ymin:badamp_conf.ymax,*]-=median(bias_diff)
                bias_adj.biasdiff-=median(bias_diff)
                n_iter=2
            endif
            if(cbad_diff gt 0) then begin
                bias_adj.biasdiff-=bias_diff
                for t=0,n_img-1 do begin
                        img_cube[badamp_conf.xmin:badamp_conf.xmax,badamp_conf.ymin:badamp_conf.ymax,t]-= bias_diff[t]
                        fl_xmin[*,t]-=bias_diff[t]
                        fl_xmin1[*,t]-=bias_diff[t]
                        fl_xmax[*,t]-=bias_diff[t]
                        fl_xmax1[*,t]-=bias_diff[t]
                        if(t lt n_img-1) then begin
                            dfl_xmin[*,t]-=(bias_diff[t+1]-bias_diff[t])
                            dfl_xmin1[*,t]-=(bias_diff[t+1]-bias_diff[t])
                            dfl_xmax[*,t]-=(bias_diff[t+1]-bias_diff[t])
                        endif
                endfor
            endif else n_iter=2
        endwhile

        if(cbad_diff gt 0 and cgood_diff ge 2) then begin
            if(keyword_set(debug)) then print,'Good readout differences, iteration 1: ',good_diff
    
            robomean,(median((dfl_xmin1-dfl_ref_min)^2,dim=1))[good_diff],3.0,0.5,dflm_xmin_good,rmi2,dfls_xmin_good
            robomean,(median((dfl_xmax-dfl_ref_max)^2,dim=1))[good_diff],3.0,0.5,dflm_xmax_good,rma2,dfls_xmax_good
            dfls_ref_min_good=stdev(median(dfl_ref_min[*,good_diff],dim=1))
            dfls_ref_max_good=stdev(median(dfl_ref_max[*,good_diff],dim=1))
            good_diff_min=fltarr(n_img-1)
            good_diff_max=fltarr(n_img-1)
            for ii=0,n_img-2 do begin
                good_diff_min[ii]=n_elements(where(abs(dfl_xmin1[*,ii]-dfl_ref_min[*,ii]) lt sqrt(2)*sqrt(dfl_ref_min[*,ii]+15.^2)/sqrt(5.0)))/1024.
                good_diff_max[ii]=n_elements(where(abs(dfl_xmax[*,ii]-dfl_ref_max[*,ii]) lt sqrt(2)*sqrt(dfl_ref_max[*,ii]+15.^2)/sqrt(5.0)))/1024.
            endfor
            gd_2 = where(good_diff_min[good_diff] gt 0.45 or good_diff_max[good_diff] gt 0.45, cgood_diff_2, compl=bd_2, ncompl=cbd_2)
;            std_thr=2.0
;            gd_2 = where(abs((median((dfl_xmin1-dfl_ref_min)^2,dim=1))[good_diff]-dflm_xmin_good) lt dfls_ref_min_good*std_thr or $
;                         abs((median((dfl_xmax-dfl_ref_max)^2,dim=1))[good_diff]-dflm_xmax_good) lt dfls_ref_max_good*std_thr, $
;                         cgood_diff_2, compl=bd_2, ncompl=cbd_2)
            good_diff_2 = good_diff[gd_2]
            if(cbd_2 gt 0) then bad_diff_2=[bad_diff,good_diff[bd_2]] else bad_diff_2=bad_diff
            cbad_diff_2=cbad_diff+cbd_2
            if(keyword_set(debug)) then print,'Good readout differences, iteration 2: ',good_diff_2
            good_readouts=[good_diff_2,good_diff_2+1L]
            good_readouts=good_readouts[sort(good_readouts)]
            ugrd=uniq(good_readouts)
            good_readouts=good_readouts[ugrd]
            cgood_readouts=n_elements(good_readouts)
    
            readouts_mask=bytarr(n_img)
            readouts_mask[good_readouts]=1
            bad_readouts=where(readouts_mask eq 0,cbad_readouts)
    
;            bias_diff=((fl_xmin[good_readouts]-fl_ref_min[good_readouts])+(fl_xmax[good_readouts]-fl_ref_max[good_readouts]))/2.0
            bias_diff=median(((fl_xmin1[*,good_readouts]-(2*fl_ref_min[*,good_readouts]-fl_ref_min1[*,good_readouts]))+$
                       (fl_xmax[*,good_readouts]-(2*fl_ref_max[*,good_readouts]-fl_ref_max1[*,good_readouts])))/2.0,dim=1)

            group_grd=[0]
            for t=1,cgood_readouts-1 do begin
                if((good_readouts[t]-good_readouts[t-1] eq 1) and (t ne cgood_readouts-1)) then $
                    group_grd=[group_grd,t] $
                else begin
                    if(t eq cgood_readouts-1) then group_grd=[group_grd,t]
                    bias_curr=mean(bias_diff[group_grd])
                    if(keyword_set(debug)) then print,'Bias adjustment in a group',group_grd,' readout ',good_readouts[t],' mean',mean(bias_diff[group_grd]),' stdev',stdev(bias_diff[group_grd])
                    ;bias_diff[group_grd]=bias_curr
                    group_grd=[t]
                endelse
            endfor
            for t=0,cgood_readouts-1 do begin
                img_cube[badamp_conf.xmin:badamp_conf.xmax,badamp_conf.ymin:badamp_conf.ymax,good_readouts[t]]-=bias_diff[t]
                bias_adj[good_readouts[t]].biasdiff-=bias_diff[t]
            endfor
            for tt=0,cbad_readouts-1 do begin
                bias_adj[bad_readouts[tt]].used=0
                img_cube[badamp_conf.xmin:badamp_conf.xmax,badamp_conf.ymin:badamp_conf.ymax,bad_readouts[tt]]=!values.f_nan
            endfor
        endif else begin
            if(cgood_diff lt 2) then begin
                message,/inf,'Less than 2 read-outs are considered acceptable, masking the bad region'
                img_cube[badamp_conf.xmin:badamp_conf.xmax,badamp_conf.ymin:badamp_conf.ymax,*]=!values.f_nan
                cbad_diff_2=n_img-1
                cgood_diff_2=0
                cbad_readouts=n_img
                cgood_readouts=0
                bad_readouts=indgen(cbad_readouts)
                good_readouts=[-1]
            endif else begin
                cbad_diff_2=0
                cgood_diff_2=n_img-1
                cbad_readouts=0
                cgood_readouts=n_img
                good_readouts=indgen(cgood_readouts)
                bad_readouts=[-1]
            endelse
        endelse

        sxaddpar,h_img,'BADAMP',badamp,' bad amplifier number (0 to 31)'
        sxaddpar,h_img,'AMP'+string(badamp,form='(i2.2)')+'NGRD',cgood_readouts,' number of good read-outs'
        sxaddpar,h_img,'AMP'+string(badamp,form='(i2.2)')+'NBRD',cbad_readouts,' number of bad read-outs'
        string_bad_readouts=(cbad_readouts eq 0)? 'none' : $
                            strjoin(string(bad_readouts+(1-keyword_set(keepfirst)),form='(i2.2)'),',')
        sxaddpar,h_img,'AMP'+string(badamp,form='(i2.2)')+'BRD',string_bad_readouts,' discarded readouts'

        if(keyword_set(verbose)) then begin
            print,'Correction attempted for the bad amplifier N'+string(badamp,form='(i2.2)')
            print,'Numbers of good and bad readouts: ',cgood_readouts,cbad_readouts,form='(a,i4,i4)'
            print,'Discarded readouts: ',string_bad_readouts
        endif
    endif else begin
        ;;;; not implemented yet for the amplifiers 8 to 15 and 24 to 31
    endelse
endif else begin
    sxaddpar,h_img,'BADAMP',-1,' no bad amplifier correction performed'
    cbad_diff_2=0
    cgood_diff_2=n_img-1
    cbad_readouts=0
    cgood_readouts=n_img
    good_readouts=indgen(cgood_readouts)
    bad_readouts=[-1]
endelse

if(keyword_set(debug)) then stop

;; the use of the IDL linear keyword seems to imply the data are
;; linear

if(not keyword_set(linear)) then begin
;; cnaw: added an environment variable so calibrations are found 2017-05-17
   calib_path = getenv('IDL_MMIRS')
   filen = calib_path+'calib_MMIRS/detector/'+detector_id+'/nonlinearity/ramp'+strcompress(string(ramp,format='(i)'),/remove)+'_nonlin.fits'
   print, 'detector_id: ', detector_id, ' ',filen, ramp
   t_cal=mrdfits(calib_path+'calib_MMIRS/detector/'+detector_id+'/nonlinearity/ramp'+strcompress(string(ramp,format='(i)'),/remove)+'_nonlin.fits',1,/silent)
    for i=0,n_img-1 do begin
        if(keyword_set(verbose)) then print,'Nonlin correction for readout N=',string(i+1,form='(i3)'),format='(a,a,%"\r",$)'
        if(n_elements(bias) gt 1) then $
            corr_img=mmirs_correct_nonlinearity(img_cube[*,*,i],calib=t_cal,/nobias) $
        else $
            corr_img=mmirs_correct_nonlinearity(img_cube[*,*,i],calib=t_cal)
        if(detector_id eq 'H2RG') then corr_img=transpose(corr_img)
        img_cube[*,*,i]=corr_img
        sat_pix=where(finite(corr_img) ne 1 and sat_mask eq 0,csat_pix)
        if(csat_pix gt 1) then sat_mask[sat_pix]=i+1
    endfor
    if(keyword_set(verbose)) then print,''
endif

if(keyword_set(debug)) then stop

sxaddpar,h_img,'NLINCORR',1-keyword_set(linear),' non-linearity correction status'

if(n_img eq 1) then begin
    img_final=img_cube[*,*,0]/ramp
endif 
if(n_img eq 2) then begin
    img_final=(img_cube[*,*,1]-img_cube[*,*,0])/ramp
endif
if(n_img gt 2) then begin
    img_final=img_cube[*,*,0]*!values.f_nan
    if(keyword_set(verbose)) then print,'Computing the difference...',format='(a,$)'
    for k=0,n_img-2 do img_cube[*,*,k]=(img_cube[*,*,k+1]-img_cube[*,*,k])/ramp
    if(keyword_set(verbose)) then print,'done'

    if(badamp ge 0) then if (cbad_diff_2 gt 0 and cgood_diff_2 ge 1) then begin ;;; filling the missing readouts for a bad amplifier with the mean values
        if(keyword_set(verbose)) then print,'Filling the bad readouts with the mean value'
        badamp_cube=img_cube[badamp_conf.xmin:badamp_conf.xmax,badamp_conf.ymin:badamp_conf.ymax,*]
        mean_bad_amp=(cgood_diff_2 eq 1)? badamp_cube[*,*,good_diff_2] : total(badamp_cube[*,*,good_diff_2],3)/cgood_diff_2
        for tt=0,cbad_diff_2-1 do img_cube[badamp_conf.xmin:badamp_conf.xmax,badamp_conf.ymin:badamp_conf.ymax,bad_diff_2[tt]]=mean_bad_amp
    endif

;        if(keyword_set(clean)) then begin
;            print,'Cleaning: '
;            nsig=10.0
;            img_med_frag=median(img_cube[*,*,0:n_img-2L],dim=3)
;            for n=0L,n_img-2L do begin
;                badpix=where(abs(img_cube[*,*,n]-img_med_frag) gt sqrt(2)*nsig*rdnoise/ramp,cbadpix)
;                if(cbadpix gt 0) then img_cube[badpix+nx*ny*n]=img_med_frag[badpix]
;                print,n,cbadpix,n_elements(img_med_frag),nsig*rdnoise/ramp
;            endfor
;        endif

;;
;; Cosmic ray removal, flag bad pixels and saturated pixels
;;
    block_size=64L
    n_bl_max=2048L/block_size
    cbadpixtot=0L
    nsig=9.0


    for iter=0,keyword_set(clean)*2 do begin
        if(keyword_set(clean) and keyword_set(verbose)) then print,'Removing cosmic ray hits. Iteration '+string(iter+1,format='(i2)')
        for n_bl=0L,n_bl_max-1L do begin
            ymin=n_bl*block_size
            ymax=ymin+block_size-1L
            idx_off=block_size*2048L*n_bl
            if(keyword_set(clean) and iter gt 0) then begin
;                print,'Cleaning: '
                img_fit_frag=img_final[*,ymin:ymax]
                for n=0L,n_img-2L do begin
                    medval=0.0; median(img_cube[*,ymin:ymax,n]-img_fit_frag)
                    stdval=robust_sigma(img_cube[*,ymin:ymax,n]-img_fit_frag)
;;;;                    badpix=where(abs(img_cube[*,ymin:ymax,n]-img_fit_frag-medval) gt sqrt((nsig*stdval)^2+(rdnoise/ramp)^2*0.),cbadpix)
                    badpix=where(abs(img_cube[*,ymin:ymax,n]-img_fit_frag-medval) gt nsig*sqrt(abs(img_fit_frag)*ramp+16.0*rdnoise^2)/ramp,cbadpix)
                    cbadpixtot=cbadpixtot+cbadpix
                    if(cbadpix gt 0) then img_cube[badpix+idx_off+nx*ny*n]=img_fit_frag[badpix]
;                    print,n_bl,n,cbadpix,n_elements(img_fit_frag),nsig*rdnoise/ramp
                endfor
            endif
            for n=0L,n_img do begin
                sat_pix=where(sat_mask[*,ymin:ymax] eq n,csat_pix)
                if(keyword_set(verbose)) then print,'block=',n_bl,' N=',n,' csat_pix=',csat_pix,format='(a,i4,a,i3,a,i11,%"\r",$)'
                if(csat_pix gt 0) then begin
                    np_max=(n eq 0)? n_img-1L : ((n-2L) > 1L)
                    idx_arr=reform((sat_pix # (lonarr(np_max)+1L))+$
                                   ((lonarr(csat_pix)+1L) # (lindgen(np_max)*nx*ny)), $
                                       csat_pix*np_max)
                    flux_arr=reform(img_cube[idx_arr+idx_off],csat_pix,np_max)
                    ;;;ierr_arr=ramp*gain/(2.0*((flux_arr*gain*ramp) > 0) + rdnoise^2)*0.0+1d ;;; something wrong here

                    if(np_max gt 2) then begin
;;;; in the read-noise limited case, the inverse weighting of diff read-outs
;;;; is described by a quadratic function as a number of readout
;;;; it is symmetric, with the maximum value of 1 in the middle of the
;;;; exposure declining on both sides. The value at 0 (and at N_max=N) is
                        w_min = 1d / (0.51550615d + n*0.24976740d)
;;;; the quadratic function is then described as:
                        a_coeff = (1d -w_min)*4d/(np_max-1d)
                        b_coeff = -a_coeff/(np_max-1d)
                        w_arr = 1.0/transpose(w_min+a_coeff*dindgen(np_max)+b_coeff*dindgen(np_max)^2)^2
                        ierr_arr = congrid(w_arr,csat_pix,np_max)
                    endif else $
                        ierr_arr = flux_arr*0+1d

                    flux_vec=((n ge 1) and (n le 3))? flux_arr : total(flux_arr*ierr_arr,2)/total(ierr_arr,2)
                    if(keyword_set(debug)) then stop
                    img_final[sat_pix+idx_off]=flux_vec[*]
                endif
            endfor
        endfor
    endfor
    if(keyword_set(verbose)) then print,''
    if(keyword_set(verbose)) then print,'Number of bad pixels: ',cbadpixtot
endif

if(detector_id eq 'H2RG' and (not keyword_set(linear))) then img_final=transpose(img_final)
if(detector_id eq 'H2RG') then img_final=refpix(img_final)
writefits,outfile,0,h_dummy
mwrfits,img_final,outfile,h_img
mwrfits,bias_adj,outfile
if(keyword_set(debug)) then stop

if(not keyword_set(compress)) then return

case compress of
    'gz': spawn,'gzip -f '+outfile
    '.gz': spawn,'gzip -f '+outfile
    'bz2': spawn,'bzip2 -f '+outfile
    'zip': spawn,'zip '+outfile+'.zip '+outfile
    ELSE: spawn,compress+' '+outfile
endcase

end

