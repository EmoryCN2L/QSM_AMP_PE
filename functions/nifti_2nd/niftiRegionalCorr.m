function [calcCorr,d,coef,fitedD] = niftiRegionalCorr(iname1, iname2, maskName, radius, oname)
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
for i = 1:numel(ind);
    [x,y,z] = ind2sub(size(im1),ind(i));
    startX = x-radius; if(startX<1) startX = 1;end;
    endX = x+radius; if(endX>sz(1)) endX = sz(1);end;
    startY = y-radius; if(startY<1) startY = 1;end;
    endY = y+radius; if(endY>sz(1)) endY = sz(2);end;
    vec1 = im1(startX:endX,startY:endY,z);
    vec2 = im2(startX:endX,startY:endY,z);
    vec1 = reshape(vec1,[numel(vec1),1]);
    vec2 = reshape(vec2,[numel(vec2),1]);
    cor = corrcoef([vec1,vec2]);
    imo(x,y,z) = cor(2,1);
end

hdro = hdr2;
hdro.fname = oname;
spm_write_vol(hdro,imo);




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


