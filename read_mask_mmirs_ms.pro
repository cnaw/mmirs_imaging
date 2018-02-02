function read_mask_mmirs_ms,filename,logfile=logfile
    openr,u,filename,/get_lun
    str=''
    count_lines=0ll
    count_hdr=0ll
    WHILE ~ EOF(u) DO BEGIN 
        readf,u,str
        count_lines=count_lines+1
        if(strmid(str,0,7) eq 'corners') then mask_corners=double((strsplit(str,/extr))[1:*])
        if(strmid(str,0,4) eq 'slit') then begin
            count_hdr=count_lines+1
            hdrline=str
        endif
    endwhile
    close,u,/all
    free_lun,u

    hdr_fields=strcompress(strsplit(hdrline,/extract),/remove_all)
    idx_slit=where(hdr_fields eq 'slit')
    idx_ra=where(hdr_fields eq 'ra')
    idx_dec=where(hdr_fields eq 'dec')
    idx_x=where(hdr_fields eq 'x')
    idx_y=where(hdr_fields eq 'y')
    idx_target=where(hdr_fields eq 'target')
    idx_object=where(hdr_fields eq 'object')
    idx_type=where(hdr_fields eq 'type')
    idx_wstart=where(hdr_fields eq 'wstart')
    idx_wend=where(hdr_fields eq 'wend')
    idx_height=where(hdr_fields eq 'height')
    idx_width=where(hdr_fields eq 'width')
    idx_offset=where(hdr_fields eq 'offset')
    idx_theta=where(hdr_fields eq 'theta')
    idx_bbox=where(hdr_fields eq 'bbox')
    idx_polygon=where(hdr_fields eq 'polygon')

    mask_data=read_asc(filename,/str,count_hdr)
    ;;mask_data=read_asc(filename,/str,44)
    ;;mask_corners=[-19.300,-33.500,19.300,33.300]
    n_slits=n_elements(mask_data[0,*])

    mask=replicate({slit:0,ra:0d,dec:0d,x:0d,y:0d,target:0,object:'none',type:'NONE',$
                    wstart:0.0,wend:0.0,height:0.0,width:0.0,offset:0.0,theta:0.0,$
                    bbox:dblarr(8),corners:mask_corners},n_slits)
    if(idx_slit   ge 0) then mask.slit=transpose(fix(mask_data[idx_slit,*]))
    if(idx_x      ge 0) then mask.x=transpose(float(mask_data[idx_x,*]))
    if(idx_y      ge 0) then mask.y=transpose(float(mask_data[idx_y,*]))
    if(idx_target ge 0) then mask.target=transpose(fix(mask_data[idx_target,*]))
    offs_col=(mask_data[7,0] eq 'BOX' or mask_data[7,0] eq 'TARGET')? 0 : 1 ;;; to solve the problem of empty object names
    if(idx_object ge 0) then mask.object=transpose(mask_data[idx_object,*])
    if(idx_type   ge 0) then mask.type=transpose(mask_data[idx_type-offs_col,*])
    if(idx_wstart ge 0) then mask.wstart=transpose(double(mask_data[idx_wstart-offs_col,*]))
    if(idx_wend   ge 0) then mask.wend=transpose(double(mask_data[idx_wend-offs_col,*]))
    if(idx_height ge 0) then mask.height=transpose(float(mask_data[idx_height-offs_col,*]))
    if(idx_width  ge 0) then mask.width=transpose(float(mask_data[idx_width-offs_col,*]))
    if(idx_offset ge 0) then mask.offset=transpose(float(mask_data[idx_offset-offs_col,*]))
    if(idx_theta  ge 0) then mask.theta=transpose(float(mask_data[idx_theta-offs_col,*]))

    if(idx_wstart eq -1 and n_elements(logfile) eq 1) then begin
        grism=def_grism(logfile)
        for i=0,n_slits-1 do begin
            wlsol=mmirs_wlsol_iguess(grism,mask[i].x)
            wstartend=poly([0,2047],wlsol)
            mask[i].wstart=wstartend[0]
            mask[i].wend=wstartend[1]
        endfor
    endif

    for i=0,n_slits-1 do begin
        mask[i].bbox=double(mask_data[idx_bbox-offs_col:idx_bbox+4-offs_col,i])
        get_coords,instr=mask_data[idx_ra[0],i]+' '+mask_data[idx_dec[0],i],crdtmp
        mask[i].ra=crdtmp[0]*15d
        mask[i].dec=crdtmp[1]
    endfor

    return,mask
end
