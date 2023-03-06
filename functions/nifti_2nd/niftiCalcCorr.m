function [calcCorr,d,coef,fitedD] = niftiCalcCorr(iname1, iname2, maskName)
%This function calculates correlation between two images masked by mask
%
%function [calcCorr,d,coef,fitedD] = niftiCalcCorr(iname1, iname2,
%maskName)
%calcCorr, correlation coeffecient between the two images
%d, a vector of raw data with two columns
%coef, coefficient of the linear fit
%fitedD, fitted Data

im1 = spm_read_vols(spm_vol(iname1));
im2 = spm_read_vols(spm_vol(iname2));


mask = spm_read_vols(spm_vol(maskName));
if(~isequal(size(im1),size(im2))||~isequal(size(im1),size(mask)))
    error('Incompatible image dimensions');
end
d = [im1(mask>0.5),im2(mask>0.5)];
coef = [ones(size(d,1),1),d(:,1)]\d(:,2);
fitedD = coef(1)+coef(2)*d(:,1);
calcCorr = corrcoef(d);
calcCorr = calcCorr(2,1);
return;


