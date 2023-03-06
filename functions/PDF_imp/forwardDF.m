function Y = forwardDF(x,mask,D)
mask_outside = 1-mask;
Y = mask.*ifftn(D.*fftn(x.*mask_outside));
return
