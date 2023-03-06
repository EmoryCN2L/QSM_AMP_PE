function niftiCalculateADC(iname,oname,bval);
%

[hdr, ft,fp,machine] = load_nii_hdr_ut(iname);
if(hdr.dime.dim(5)~=2)
    error('Should be a 4D volume containing both b0 and higher b image');
end;

imgB0 = niftiScaleImage(load_nii_img(hdr,ft,fp,machine,1),hdr);
imgB2 = niftiScaleImage(load_nii_img(hdr,ft,fp,machine,2),hdr);

mask = (imgB0>0)&(imgB2>0);

imgOutMasked = log(imgB2(mask)./imgB0(mask))/-bval;

ohdr = hdr;
ohdr.dime.datatype = niftiType('float32');
ohdr.dime.bitpix = 32;
ohdr.dime.dim(1) = 4;
ohdr.dime.dim(5) = 1;

ohdr.dime.scl_slope = 1;
ohdr.dime.scl_interp = 0;

ohdr.hist.magic = 'n+1';

fid = fopen(oname,'wb');
save_nii_hdr(ohdr,fid);
niftiPadHeader(ohdr,fid);


imgOut = zeros(size(imgB0));
imgOut(mask) = imgOutMasked;
write_nii_segment(fid,imgOut,ohdr);
fclose(fid);

