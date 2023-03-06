function niftiVoxelScale(iname,outname,factor);
%niftiVoxelScale(iname,outname,factor);
%scale the image dimension

[hdr, ft, fp, machine] = load_nii_hdr_ut(iname);
img = load_nii_img(hdr,ft,fp,machine);
hdr.dime.pixdim(2:4) = hdr.dime.pixdim(2:4)*factor;
niftiSaveNii(hdr,img,outname,ft);

return;
