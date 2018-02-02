function get_mmirs_period,filename,exten=exten,$
    year=year,month=month,suffix=suffix

;;if(n_elements(exten) ne 1) then exten=1
;;
;; The code above crashes with 2010 data because
;; exten is undefined and (n_elements(exten) = 1)
;;
;; cnaw 2017-07-14
;;
fits_open,filename, fcb
extend = 0
if(fcb.nextend > 0) then extend = 1


h=headfits(filename,exten=exten,/silent)
dateobs=sxpar(h,'DATE-OBS',count=count)
if(count eq 0) then begin
    year=2010
    month=0
    suffix='A'
endif else begin
    year=float(strmid(dateobs,0,4))
    month=float(strmid(dateobs,5,2))
    suffix=(month le 6)? 'A' : 'B'
endelse

period = string(year,format='(i4.4)')+suffix
print,' get_mmirs_period: period :',period
;;
;; Added this because of "too many files open"
fits_close,fcb
return,period

end
