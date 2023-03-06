function [calcCorr,d,coef,fitedD] = niftiFitResidual(iname1, iname2, maskName, oname)
%This function calculates local correlation between two images in a region specified
%function niftiRegionalCorr(iname1, iname2,maskName, oname)
%For each voxel, voxels within a specified radius will be used for the calculation of correlation coeffecient
%iname1, first input image
%iname2, second input image
%maskName, mask specifying where to calculate the regional correlation
%radius, radius of voxels used for the calculation
%oname, outputname

hdr1 = spm_vol(iname1);
hdr2 = spm_vol(iname2);
hdrMask = spm_vol(maskName);

if(~isCompatible(hdr1,hdr2)||~isCompatible(hdr1,hdrMask));
    error('Incompatible images');
end;

im1 = spm_read_vols(hdr1);
im2 = spm_read_vols(hdr2);
imMask = spm_read_vols(hdrMask);

ind = find(imMask>0.5);

sz = size(im1);
imo = im1*0;

[calcCorr,d,coef,fitedD] = niftiCalcCorr(iname1, iname2, maskName);

imo(ind) = abs(im2(ind)-(im1(ind)*coef(2)+coef(1)));


hdro = hdr2;
hdro.fname = oname;
spm_write_vol(hdro,imo);

return;



return;

function r = isCompatible(hdr1,hdr2)
if(~isequal(hdr1.dim,hdr2.dim))
    r = false;
    return;
end;

if(any(abs(hdr1.mat-hdr2.mat)>0.1))
    r = false;
    return;
end;
r = true;

return;


