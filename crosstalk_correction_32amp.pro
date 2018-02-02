function crosstalk_correction_32amp,img,sat_val=sat_val, h_img=h_img

corr_img=fltarr(2048,2048)
s_img=size(img)
if(s_img[1] ne 2048 or s_img[2] ne 2048) then begin
    message,/inf,'Only 2048x2048 images are supported. Returning original image'
    return,corr_img
endif


mask_img=bytarr(1024,1024)
wght=[0.0,0.6,0.8,1.0,1.0,1.0,1.0,1.0]

corr_val=-25.0 ;; default

;; cnaw 2017-06-07
;;corr_val = -37.5 ;;
;;corr_val = -50.0 ;; makes a difference
;;corr_val = +50.0 ;; makes no difference
if(n_elements(sat_val) ne 1) then sat_val=25000.0

if(keyword_set(h_img)) then begin
   sxaddpar, h_img,'XTALKVAL',corr_val,' cross-talk correction'
   sxaddpar, h_img,'XTALKSAT',sat_val,' cross-talk saturation level'
endif

;;; Q1
mask_img=mask_img*0b
img_q=img[1024:*,1024:*]
corr_img_q=fltarr(1024,1024)
high_img=where(img_q ge sat_val,chigh_img)
if(chigh_img gt 0) then begin
    mask_img[high_img]=1
;;    mask_img=byte(smooth(mask_img*20.0,0))
    for i=1,7 do begin
        mask_img1=corr_img_q*0.0
        mask_img1[where(mask_img gt 0)]+=corr_val*wght[i]
        mask_img1=shift(mask_img1,128*i,0)
        mask_img1[0:128*i-1,*]=0.0
        corr_img_q+=mask_img1
        mask_img1=corr_img_q*0.0
        mask_img1[where(mask_img ne 0)]+=corr_val*wght[i]
        mask_img1=shift(mask_img1,-128*i,0)
        mask_img1[1023-128*i:*,*]=0
        corr_img_q+=mask_img1
    endfor
    corr_img[1024:*,1024:*]=corr_img_q
endif

;;; Q2
mask_img=mask_img*0b
img_q=img[0:1023,1024:*]
corr_img_q=fltarr(1024,1024)
high_img=where(img_q ge sat_val,chigh_img)
if(chigh_img gt 0) then begin
    mask_img[high_img]=1
;    for i=1,7 do corr_img_q[where(shift(mask_img,0,128*i) ne 0)]+=corr_val
    for i=1,7 do begin
        ds=(i eq 1)? 2 : 0
        mask_img1=corr_img_q*0.0
        mask_img1[where(mask_img gt 0)]+=corr_val*wght[i]
        mask_img1=shift(mask_img1,0,128*i+ds)
        mask_img1[*,0:((128*i-1+ds)<1023)]=0.0
        corr_img_q+=mask_img1
        mask_img1=corr_img_q*0.0
        mask_img1[where(mask_img ne 0)]+=corr_val*wght[i]
        mask_img1=shift(mask_img1,0,-128*i+ds)
        mask_img1[*,((1023-128*i+ds)>0):*]=0
        corr_img_q+=mask_img1
    endfor
    corr_img[0:1023,1024:*]=corr_img_q
endif

;;; Q3
mask_img=mask_img*0b
img_q=img[0:1023,0:1023]
corr_img_q=fltarr(1024,1024)
high_img=where(img_q ge sat_val,chigh_img)
if(chigh_img gt 0) then begin
    mask_img[high_img]=1
    for i=1,7 do begin
        ds=(i eq 1)? -2 : 0
        mask_img1=corr_img_q*0.0
        mask_img1[where(mask_img ne 0)]+=corr_val*wght[i]
        mask_img1=shift(mask_img1,128*i+ds,0)
        mask_img1[0:((128*i-1+ds)<1023),*]=0.0
        corr_img_q+=mask_img1
        mask_img1=corr_img_q*0.0
        mask_img1[where(mask_img ne 0)]+=corr_val*wght[i]
        mask_img1=shift(mask_img1,-128*i+ds,0)
        mask_img1[((1023-128*i+ds)>0):*,*]=0
        corr_img_q+=mask_img1
    endfor
    corr_img[0:1023,0:1023]=corr_img_q
endif

;;; Q4
mask_img=mask_img*0b
img_q=img[1024:*,0:1023]
corr_img_q=fltarr(1024,1024)
high_img=where(img_q ge sat_val,chigh_img)
if(chigh_img gt 0) then begin
    mask_img[high_img]=1
;    for i=1,7 do corr_img_q[where(shift(mask_img,0,128*i) ne 0)]+=corr_val
    for i=1,7 do begin
        ds=(i eq 1)? 0 : 0
        mask_img1=corr_img_q*0.0
        mask_img1[where(mask_img gt 0)]+=corr_val*wght[i]
        mask_img1=shift(mask_img1,0,128*i+ds)
        mask_img1[*,0:((128*i-1+ds)<1023)]=0.0
        corr_img_q+=mask_img1
        mask_img1=corr_img_q*0.0
        mask_img1[where(mask_img ne 0)]+=corr_val*wght[i]
        mask_img1=shift(mask_img1,0,-128*i+ds)
        mask_img1[*,((1023-128*i+ds)>0):*]=0
        corr_img_q+=mask_img1
    endfor
    corr_img[1024:*,0:1023]=corr_img_q
endif

return,corr_img

end
