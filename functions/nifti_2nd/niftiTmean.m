function niftiTmean(iname,oname);
%compute the mean of images along the time dimension

[ihdr,ift,ifp,imachine] = load_nii_hdr_ut(iname);

ohdr = ihdr;
ohdr.hist.magic = 'n+1';%force the output to be .nii
ohdr.dime.vox_offset = 352;

ohdr.dime.dim(1) = 4;
ohdr.dime.dim(5) = 1;%mean image
ohdr.dime.scl_slope = 1;
ohdr.dime.scl_inter = 0;
ohdr.datatype = 16;%single float
ohdr.bitpix = 32;

oImg = zeros(ihdr.dime.dim(2:4));

for i = 1:ihdr.dime.dim(5);
  img = load_nii_img(ihdr,ift,ifp,imachine,i);
  oImg = oImg+niftiScaleImage(img,ihdr)/ihdr.dime.dim(5);
  disp([num2str(i),'completed']);
end;

ofid = fopen(oname,'wb');
save_nii_hdr(ohdr,ofid);
niftiPadHeader(ohdr,ofid);;

write_nii_segment(ofid,oImg,ohdr);

fclose(ofid);
return;

 