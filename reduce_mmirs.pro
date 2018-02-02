function get_filename,str,suffix=suffix
   if(n_elements(suffix) ne 1) then suffix=''
   return,[strsplit(str,',',/extract)+suffix]
end


pro reduce_mmirs,logfile,verbose=verbose

suffix='.fix.fits'

log = readlog(logfile)
print,'logfile ',logfile
rawdir = sxpar(log,'RAW_DIR')
rdir = sxpar(log,'R_DIR')
wdir = sxpar(log,'W_DIR')
;file_mkdir,wdir
spawn,'mkdir -p '+rdir
spawn,'mkdir -p '+wdir

grism = sxpar(log,'GRISM')
filter = sxpar(log,'FILTER')
slit_id = def_slit(logfile)

rawext = sxpar(log,'RAWEXT',count=cntre)
if(cntre eq 1) then suffix=suffix+rawext else rawext=''

sci_pref = sxpar(log,'SCI')
dithpos = sxpar(log,'DITHPOS')
dark_sci_all = sxpar(log,'DARKSCI')

sci_pref2 = sxpar(log,'SCI2',count=cntval)
diffmode = (cntval eq 1)? 1 : 0

arc_all = sxpar(log,'ARC')
dark_arc_all = sxpar(log,'DARKARC')

flat_all = sxpar(log,'FLAT')
dark_flat_all = sxpar(log,'DARKFLAT')

proc_stages=bytarr(9)
for i=0,n_elements(proc_stages)-1 do proc_stages[i]=sxpar(log,'S'+string(i+1,format='(i2.2)')+'PROC')

if(proc_stages[0] eq 1) then begin ;;;; dark subtraction and mask file copying
    conf_str = {rawdir:rawdir,$
                rdir:rdir,$
                wdir:wdir,$
                sci:get_filename(sci_pref),$
                sci_dark:get_filename(dark_sci_all),$
                arc:get_filename(arc_all),$
                arc_dark:get_filename(dark_arc_all),$
                flat:get_filename(flat_all),$
                flat_dark:get_filename(dark_flat_all),$
                misc:[''],misc_dark:get_filename(dark_sci_all)}

    if(diffmode eq 1) then begin
        conf_str.misc=[sci_pref2]
    endif

    status = check_input_files_mmirs(conf_str,suffix=suffix)
    print,'reduce_mmirs: do preprocessed images exist ? status ', status,' ', suffix
    if(((status gt 1) and (strlen(conf_str.misc[0]) gt 0)) or $
       ((status mod 16) gt 1)) then begin
        message,/inf,'Some files require pre-processing. Trying to run mmfixen_nonlin.'
        status_raw = check_input_files_mmirs(conf_str,/raw)
        print, 'reduce_mmirs: status_raw:', status_raw,(status_raw mod 16) 
        if((status_raw mod 16) le 1) then begin
            preproc_mmirs,conf_str,rawext=rawext,/verbose

            if(((status_raw mod 2) eq 0) and (slit_id eq 'mos')) then begin
                file_delete,(rdir+get_filename(sci_pref,suffix='.msk'))[0],/allow_nonexistent
                file_copy,(rawdir+get_filename(sci_pref,suffix='.msk'))[0],(rdir+get_filename(sci_pref,suffix='.msk'))[0],/overwrite
            endif

            status = check_input_files_mmirs(conf_str,suffix=suffix)
        endif else begin
            message,/inf,'Some files are missing in the raw archive. Cannot continue. Status='+string(status,format='(i)')
            return
        endelse
    endif

    if(slit_id eq 'mos') then begin
        if((status mod 2) eq 1) then begin
            message,/inf,'Mask definition file for '+conf_str.sci+' is missing. Cannot continue'
            return
        endif
        file_delete,wdir+'mask_mos.txt',/allow_nonexistent
        file_copy,(rdir+get_filename(sci_pref,suffix='.msk'))[0],wdir+'mask_mos.txt',/overwrite
    endif else begin
        file_delete,wdir+'mask_mos.txt',/allow_nonexistent
        file_copy,'calib_MMIRS/LS/'+slit_id+'.txt',wdir+'mask_mos.txt',/overwrite
    endelse

    print,'Subtracting dark frames'
    subtract_dark_mmirs,conf_str,suffix=suffix,/median
    print,'done - status OK'
endif

if(proc_stages[1] eq 1) then begin
    print,'Creating slit distortion map'
    distortion_mmirs_ms,logfile
    print,'done - status OK'
endif

if(proc_stages[2] eq 1) then begin
    print,'Creating and applying normalized flat field'
    n_sl = sxpar(log,'FFNSLIT',count=n_n_sl)
    norm_slit=((n_sl eq -2) or (n_n_sl eq 0))? [0,0] : [n_sl]
    flat_fielding_mmirs_ms,logfile,['obj'],norm_slit=norm_slit,$
        sub_sc_flat=sxpar(log,'FFSCFLAT'),sub_sc_sci=sxpar(log,'FFSCSCI')
    flat_fielding_mmirs_ms,logfile,['arc','flat'],dymask=0,norm_slit=norm_slit,$
        sub_sc_flat=sxpar(log,'FFSCFLAT'),sub_sc_sci=sxpar(log,'FFSCSCI')
    if(diffmode eq 1) then begin
        flat_fielding_mmirs_ms,logfile,sci_pref2,norm_slit=norm_slit,$
            sub_sc_flat=sxpar(log,'FFSCFLAT'),sub_sc_sci=sxpar(log,'FFSCSCI')
        sci1=readfits(wdir+'obj_ff.fits',hsci1)
        sci2=readfits(wdir+sci_pref2+'_ff.fits',hsci2)
        writefits,wdir+'obj_diff_ff.fits',sci1-sci2,hsci1
    endif

    print,'done - status OK'
    print,'Extracting 2D slits'
    extract_2d_slits,logfile,'obj' 
    extract_2d_slits,logfile,/nflat
    extract_2d_slits,logfile,'arc'
    print,'done - status OK'

    if(diffmode eq 1) then begin
        print,'Processing difference image'
        flat_fielding_mmirs_ms,logfile,sci_pref2,norm_slit=norm_slit,$
            sub_sc_flat=sxpar(log,'FFSCFLAT'),sub_sc_sci=sxpar(log,'FFSCSCI')
        sci1=readfits(wdir+'obj_ff.fits',hsci1)
        sci2=readfits(wdir+sci_pref2+'_ff.fits',hsci2)
        writefits,wdir+'obj_diff_ff.fits',sci1-sci2,hsci1
        extract_2d_slits,logfile,'obj_diff'
        print,'done - status OK'
    endif


endif

if(proc_stages[3] eq 1) then begin
    create_wavesol_ms,logfile,ndeg=sxpar(log,'WLNDEG'),y_ndeg=sxpar(log,'WLYNDEG'),$
        plot=sxpar(log,'WLPLOT'),$
        debug=sxpar(log,'WLDEBUG'),oh=fix(sxpar(log,'WLOH'))
    if(sxpar(log,'WLADJ') eq 1) then begin
        stages_adj_wl = (slit_id eq 'mos')? [1,1] : [0,1]
        adjust_wavesol_ms,logfile,stages=stages_adj_wl,$
            debug=sxpar(log,'WLDEBUG')
    endif
endif

if(proc_stages[4] eq 1) then begin ;;; sky subtraction
    hsci1=headfits(wdir+'obj_ff.fits')

    ssdim = (slit_id eq 'mos')? sxpar(log,'SSDIM') : 2

    ssadjwl = sxpar(log,'SSADJWL')
    ;;; everyn=100 makes no sense, increased to 1500
    create_sky_mmirs_ms,logfile,'obj',adj=ssadjwl,dim=ssdim,everyn=500,debug=sxpar(log,'SSDEBUG'),nomask=(sxpar(hsci1,'AMP20NBRD') eq 0)? 1 : 0
    sub_sky_mmirs_ms,logfile,'obj',adj=ssadjwl

    if(diffmode eq 1) then begin
        hsci2=headfits(wdir+sci_pref2+'_ff.fits')
        create_sky_mmirs_ms,logfile,'obj_diff',adj=ssadjwl,dim=ssdim,everyn=500,debug=sxpar(log,'SSDEBUG'),nomask=((sxpar(hsci1,'AMP20NBRD') eq 0) and (sxpar(hsci2,'AMP20NBRD') eq 0))? 1 : 0
        sub_sky_mmirs_ms,logfile,'obj_diff',adj=ssadjwl
    endif
    analyse_sky,logfile,bin_s=0.1
    analyse_sky,logfile,bin_s=0.1,/raw
endif

if(proc_stages[5] eq 1) then begin ;;; linearisation
    print,'Performing linearisation of 2D spectra'
    linearisation,logfile,'obj',adj=sxpar(log,'LINADJWL') ;;;,$ ;;;; non sky subtracted
;;;        subskybox=sxpar(log,'LINSSBOX'),subskytarget=sxpar(log,'LINSSTRG')
    linearisation,logfile,'obj-sky',adj=sxpar(log,'LINADJWL'),$
        subskybox=sxpar(log,'LINSSBOX'),subskytarget=sxpar(log,'LINSSTRG')
    linearisation,logfile,'arc',adj=sxpar(log,'LINADJWL')
    if(diffmode eq 1) then begin
        crhp_reject_ms,logfile
        linearisation,logfile,'obj_diff',adj=sxpar(log,'LINADJWL'),$
            subskybox=sxpar(log,'LINSSBOX'),subskytarget=sxpar(log,'LINSSTRG')
        linearisation,logfile,'obj_diff-sky',adj=sxpar(log,'LINADJWL'),$
            subskybox=sxpar(log,'LINSSBOX'),subskytarget=sxpar(log,'LINSSTRG'),/usebadpixmask
    endif
endif

if(proc_stages[6] eq 1) then begin
    print,'Extracting 1D spectra'
    extract_1d_slits,logfile,'obj',dith_from_box=sxpar(log,'EXTDFBOX'),box_exp=sxpar(log,'EXTBOXEX')
    extract_1d_slits,logfile,'obj-sky',dith_from_box=sxpar(log,'EXTDFBOX'),box_exp=sxpar(log,'EXTBOXEX')
    if(diffmode eq 1) then begin
        extract_1d_slits,logfile,'obj_diff',dith_from_box=sxpar(log,'EXTDFBOX'),box_exp=sxpar(log,'EXTBOXEX'),/diffmode
        extract_1d_slits,logfile,'obj_diff-sky',dith_from_box=sxpar(log,'EXTDFBOX'),box_exp=sxpar(log,'EXTBOXEX'),/diffmode
    endif
endif


;; process telluric stars up to extraction

if(proc_stages[7] eq 1) then begin
    n_tel=0
    ttcnt=1
    while(ttcnt gt 0) do begin
        stmp=sxpar(log,'STAR'+string(n_tel+1,format='(i2.2)'),count=ttcnt)
        if(ttcnt eq 1) then n_tel=n_tel+1        
    endwhile
    for i=1,n_tel do begin
        tel_list=get_telluric_list(logfile,n_telluric=i)
        dark_tel_list=sxpar(log,'DARKST'+string(i,format='(i2.2)'))
        tel_type=sxpar(log,'STTYPE'+string(i,format='(i2.2)'),count=ttcnt)
        if(ttcnt ne 1) then tel_type='a0v'
        conf_str_tel = {rawdir:rawdir,$
                        rdir:rdir,$
                        wdir:wdir,$
                        sci:'',sci_dark:[''],$
                        arc:'',arc_dark:[''],$
                        flat:'',flat_dark:[''],$
                        misc:tel_list,$
                        misc_dark:get_filename(dark_tel_list)}

        status = check_input_files_mmirs(conf_str_tel,suffix=suffix)
        if(status ge 16) then begin
            message,/inf,'Some telluric files require pre-processing. Trying to run mmfixen_nonlin.'
            status_raw = check_input_files_mmirs(conf_str_tel,/raw)
            if(status_raw lt 16) then begin
                preproc_mmirs,conf_str_tel,rawext=rawext,/verbose
            endif else begin
                message,/inf,'Some telluric files are missing in the raw archive. Cannot continue'
                return
            endelse
        endif

        print,'Subtracting dark frames from telluric spectra, STAR'+string(i,format='(i2.2)')
        subtract_dark_mmirs,conf_str_tel,suffix=suffix,/median
        print,'done - status OK'

        print,'Creating and applying normalized flat field for telluric spectra'
        n_sl = sxpar(log,'FFNSLIT',count=n_n_sl)
        norm_slit=((n_sl eq -1) or (n_n_sl eq 0))? [0,0] : [n_sl]
        flat_fielding_mmirs_ms,logfile,tel_list,norm_slit=norm_slit,$
            sub_sc_flat=sxpar(log,'FFSCFLAT'),sub_sc_sci=sxpar(log,'FFSCSCI')
        print,'done - status OK'
        print,'Identifying slits with telluric spectra:'
        tel_slit_id=telluric_slit_id(logfile,n_telluric=i)
        print,'done - status OK'
        print,'Extracting 2D slits from telluric spectra'

        mask=read_mask_mmirs_ms(wdir+'mask_mos.txt')
        n_slits=n_elements(mask)

        tel_slit_files=strarr(n_slits)+tel_list[0]
        for j=0,n_elements(tel_slit_id)-1 do $
            if(tel_slit_id[j].slit gt 0) then $
                tel_slit_files[tel_slit_id[j].slit-1]=tel_slit_id[j].obs_id
        extract_2d_slits,logfile,tel_slit_files,$
            out_img_type='star_tel_'+string(i,format='(i2.2)')
        print,'done - status OK'

        linearisation,logfile,'star_tel_'+string(i,format='(i2.2)'),$
            adj=sxpar(log,'LINADJWL'),$
            subskybox=sxpar(log,'LINSSBOX'),subskytarget=sxpar(log,'LINSSTRG'),/telluric
        linearisation,logfile,['flatn'],adj=sxpar(log,'LINADJWL')

        extract_1d_slits,logfile,'star_tel_'+string(i,format='(i2.2)'),/detect,/optimal
        define_telluric_correction,logfile,n_tel=i,modelstar=tel_type ;,vr=vr
    endfor
endif

;; Correct object spectra for telluric absorption

if(proc_stages[8] eq 1) then begin
    n_tel=0
    ttcnt=1
    while(ttcnt gt 0) do begin
        stmp=sxpar(log,'STAR'+string(n_tel+1,format='(i2.2)'),count=ttcnt)
        if(ttcnt eq 1) then n_tel=n_tel+1        
    endwhile
    noextr=(slit_id eq 'mos')? 0 : 1
    print,'reduce_mmirs noextr, slit_id : ', noextr, ' ',slit_id
    correct_telluric,logfile,range(1,n_tel),diffmode=diffmode,noextr=noextr
endif

end
