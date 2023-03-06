function niftiNoneZeroMask(iname,oname);
%

[ihdr,ift,ifp,imachine] = load_nii_hdr_ut(iname);


ohdr = ihdr;

ohdr.dime.dim(1) = 4;
ohdr.dime.dim(5) = 1;
ohdr.dime.vox_offset = 352;
ohdr.dime.scl_slope = 1;
ohdr.dime.scl_inter = 0;
ohdr.dime.datatype = niftiType('uint8');
ohdr.dime.bitpix = 8;

oImg  = ones(ihdr.dime.dim(2:4));

for i = 1:ihdr.dime.dim(5);
  img = load_nii_img(ihdr,ift,ifp,imachine,i);
  oImg(img==0)=0;
  disp([num2str(i),'completed']);
end;


fid = fopen(oname,'wb');
save_nii_hdr(ohdr,fid);
niftiPadHeader(ohdr,fid);
write_nii_segment(fid,oImg,ohdr);
fclose(fid);

return

  