function check_input_files_mmirs,conf_str,res_str=res_str,suffix=suffix,raw=raw

if(n_elements(suffix) ne 1) then suffix='.fix.fits'

res_str=conf_str
status=0

tstdir = (keyword_set(raw))? conf_str.rawdir : conf_str.rdir

if(keyword_set(raw)) then suffix='.fits*'

if(~file_test(tstdir+conf_str.sci+'.msk')) then status=status+1

if((strlen(conf_str.sci[0]) gt 0) and (strlen(conf_str.sci_dark[0]) gt 0)) then begin
    ff = 0
    res_str.sci = (file_test(tstdir+conf_str.sci[0]+suffix))? '1' : '0'
    ff+= (file_test(tstdir+conf_str.sci+suffix))? 0 : 1
    for i=0,n_elements(conf_str.sci_dark)-1 do begin
        res_str.sci_dark[i] = file_test(tstdir+conf_str.sci_dark[i]+suffix)? '1' : '0'
;;        print,'check_input_files_mmirs: ',tstdir+conf_str.sci_dark[i]+suffix
    ff+= file_test(tstdir+conf_str.sci_dark[i]+suffix)? 0 : 1
    endfor
    if(ff ne 0) then status=status+2
;;    print, 'check_input_files_mmirs: sci_dark : ff, status', ff, status
endif else status=status+2

if((strlen(conf_str.arc[0]) gt 0) and (strlen(conf_str.arc_dark[0]) gt 0)) then begin
    ff = 0
    for i=0,n_elements(conf_str.arc)-1 do begin
        res_str.arc[i] = (file_test(tstdir+conf_str.arc[i]+suffix))? '1' : '0'
;;        print,'check_input_files_mmirs: ',tstdir+conf_str.arc[i]+suffix
        ff+= (file_test(tstdir+conf_str.arc[i]+suffix))? 0 : 1
    endfor
    for i=0,n_elements(conf_str.arc_dark)-1 do begin
        res_str.arc_dark[i] = file_test(tstdir+conf_str.arc_dark[i]+suffix)? '1' : '0'
;;        print,'check_input_files_mmirs: ',tstdir+conf_str.arc_dark[i]+suffix
        ff+= file_test(tstdir+conf_str.arc_dark[i]+suffix)? 0 : 1
    endfor
    if(ff ne 0) then status=status+4
;;    print, 'check_input_files_mmirs: arc ', ff,status
endif else status=status+4

if((strlen(conf_str.flat[0]) gt 0) and (strlen(conf_str.flat_dark[0]) gt 0)) then begin
    ff = 0
    for i=0,n_elements(conf_str.flat)-1 do begin
        res_str.flat[i] = (file_test(tstdir+conf_str.flat[i]+suffix))? '1' : '0'
;;        print,'check_input_files_mmirs: ',tstdir+conf_str.flat[i]+suffix
        ff+= (file_test(tstdir+conf_str.flat[i]+suffix))? 0 : 1
    endfor
    for i=0,n_elements(conf_str.flat_dark)-1 do begin
        res_str.flat_dark[i] = file_test(tstdir+conf_str.flat_dark[i]+suffix)? '1' : '0'
;;        print,'check_input_files_mmirs: ',tstdir+conf_str.flat_dark[i]+suffix
        ff+= file_test(tstdir+conf_str.flat_dark[i]+suffix)? 0 : 1
    endfor
    if(ff ne 0) then status=status+8
;;    print, 'check_input_files_mmirs: flat ', ff,status
endif else status=status+8

if((strlen(conf_str.misc[0]) gt 0) and (strlen(conf_str.misc_dark[0]) gt 0)) then begin
    ff = 0
    for i=0,n_elements(conf_str.misc)-1 do begin
        res_str.misc[i] = (file_test(tstdir+conf_str.misc[i]+suffix))? '1' : '0'
;;        print,'check_input_files_mmirs: ',tstdir+conf_str.misc[i]+suffix
        ff+= (file_test(tstdir+conf_str.misc[i]+suffix))? 0 : 1
    endfor
    for i=0,n_elements(conf_str.misc_dark)-1 do begin
        res_str.misc_dark[i] = file_test(tstdir+conf_str.misc_dark[i]+suffix)? '1' : '0'
;;        print,'check_input_files_mmirs: ',tstdir+conf_str.misc_dark[i]+suffix
        ff+= file_test(tstdir+conf_str.misc_dark[i]+suffix)? 0 : 1
    endfor
    if(ff ne 0) then status=status+16
;;    print, 'check_input_files_mmirs:misc, misc_dark ', ff, status
    
endif else status=status+16


return,status

end
