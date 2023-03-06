function niftiCrossMerge(iname1,iname2,oname,lowstackslices)
%This function cross merges two volumes by taking the firt lowstackslices
%from the iname1 and other higher slices from iname2
%
%function niftiCrossMerge(iname1,iname2,oname,lowstackslices)
%
%


hdr1 = spm_vol(iname1);
hdr2 = spm_vol(iname2);

if(~isequal(hdr1.mat,hdr2.mat))
    error('Incompatible geometry');
end;

img1 = spm_read_vols(hdr1);
img2 = spm_read_vols(hdr2);

oimg = img1;
oimg(:,:,(lowstackslices+1):end) = img2(:,:,(lowstackslices+1):end);

ohdr = hdr1;
ohdr.fname = oname;

spm_write_vol(ohdr,oimg);
return;
