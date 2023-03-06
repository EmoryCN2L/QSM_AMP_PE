function x = backwardDF(Y,mask,D)
mask_outside = 1-mask;
x = mask_outside.*ifftn(conj(D).*fftn(Y.*mask));
return
